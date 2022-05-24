#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <tuple>
#include <random>
#include <algorithm>
#include <gurobi_c++.h>
#include <lemon/list_graph.h>
#include <lemon/kruskal.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>


using namespace std;
using namespace lemon;

// Defining types
typedef lemon::ListGraph Graph;
typedef lemon::ListGraph::Node Node;
typedef lemon::ListGraph::Edge Edge;
typedef lemon::ListGraph::EdgeMap<double> EdgeMapDouble;
typedef lemon::ListGraph::EdgeMap<bool>   EdgeMapBool;
typedef lemon::FilterEdges<Graph,EdgeMapBool> Tree;
typedef lemon::FilterEdges<Graph,EdgeMapBool>::NodeMap<bool>    TreeNodeMapBool;
typedef lemon::FilterNodes<Tree, vector<Node>> SubTree;
typedef tuple<int, int> uv;
typedef tuple<int, int, int> ouv;


auto seed = 0xf0da5e;

// Defining Graph elements
Graph graph;
int n, m; // Number of nodes and edges, respectively
long unsigned int k;
Node root;
vector<Node> nodes;
vector<Edge> edges;
EdgeMapDouble lengths(graph);
vector<vector<double>> requirements = {};

// Defining Solution elements
EdgeMapBool edges_tree(graph);
Tree tree(graph, edges_tree);
TreeNodeMapBool in_some_cluster(tree, false);
vector<vector<Node>> clusters;
vector<Node>* current_cluster = nullptr;

// Solver Elements
GRBEnv env(false);


void readInstance(char filename[])
{
    // Stream and temporary vars
    ifstream instanceFile(filename);
    string e1, e2, e3;
    uint16_t u, v;
    double value;

    // Get graph size (n for nodes and m for edges) and resizing everything
    instanceFile >> e1 >> e2;
    n = stoi(e1);
    m = stoi(e2);
    graph.reserveNode(n);
    nodes.resize(n);
    requirements.resize(n, vector<double>(n));
    graph.reserveEdge(m);
    edges.resize(m);

    // Adding nodes
    for (auto i=0; i<n; i++)
    {
        nodes[i] = graph.addNode();
    }

    // Adding edges and length
    for (auto i=0; i<m; i++)
    {
        instanceFile >> e1 >> e2 >> e3;
        u = stoi(e1);
        v = stoi(e2);
        value = stod(e3);
        edges[i] = graph.addEdge(nodes[u], nodes[v]);
        lengths[edges[i]] = value;
    }
    // Adding requirements
    for (auto i=0; i<n; i++)
    {
        requirements[i][i] = 0;
        for (auto j=i+1; j<n; j++)
        {
            instanceFile >> e1;
            value = stod(e1);
            requirements[i][j] = value;
            requirements[j][i] = value;
        }
    }
    instanceFile.close();
}

void printEdgesTree()
{
    for (auto edge : edges)
    {
        if (edges_tree[edge])
        {
            cout << graph.id(graph.u(edge)) << " " << graph.id(graph.v(edge)) << endl;
        }
    }
}

double generateInitialSolution()
{
    return kruskal(graph, lengths, edges_tree);
}

double calculateObjective()
{
    double cost = 0;
    for (auto u : nodes)
    {
        Dijkstra<Tree, EdgeMapDouble> dij(tree, lengths);
        dij.init();
        dij.addSource(u);
        // Avoiding processing the u-u edge
        dij.processNextNode();
        int counter = 0;
        while (!dij.emptyQueue())
        {
            counter++;
            Node v = dij.processNextNode();
            // cout << graph.id(u) << " " << graph.id(v) << endl;
            double distance = dij.dist(v);
            double requirement = requirements[tree.id(u)][tree.id(v)];
            cost += distance * requirement;
        }
        // cout << "esges visited: " << counter << endl;
    }
    return cost/2;
}

void printClusters()
{
    for (auto cluster : clusters)
    {
        cout << "{ ";
        for (auto node : cluster)
        {
            cout << graph.id(node) << " ";
        }
        cout << "}" << endl;
    }
}

void addToCluster(Node n)
{
    if (current_cluster == nullptr)
    {
        current_cluster = new vector<Node>();
    }
    in_some_cluster[n] = true;
    current_cluster->push_back(n);
    if (current_cluster->size() == k)
    {
        clusters.push_back(*current_cluster);
        // printClusters();
        current_cluster = nullptr;
    }
}

bool divideTree(Node n, Node up)
{
    if (n==INVALID) return false;

    if (!in_some_cluster[n])
    {
        // cout << "Node " << tree.id(n) << " added for the first time" << endl;
        addToCluster(n);
    }

    for (Tree::EdgeIt it(tree); it != INVALID; ++it)
    {
        Node next_node = tree.oppositeNode(n, it);
        if (next_node == up) continue;
        bool node_exists = divideTree(next_node, n);
        if (node_exists)
        {
            if (current_cluster != nullptr &&
                find(current_cluster->begin(), current_cluster->end(), n) == current_cluster->end())
            {
                // cout << "Node " << tree.id(n) << " added for connection" << endl;
                addToCluster(n);
            }
        }
    }
    if (up == INVALID && current_cluster != nullptr)
    {
        clusters.push_back(*current_cluster);
        current_cluster = nullptr;
    }
    return true;
}



void printTree(Node n, Node up)
{
    for (Tree::IncEdgeIt e(tree, n); e != INVALID; ++e)
    {
        if (tree.oppositeNode(n,e) != up)
        {
            cout << tree.id(n) << " " << tree.id(tree.oppositeNode(n,e)) << endl;
            printTree(tree.oppositeNode(n,e), n);
        }
    }
}


void addRequirements(int base_index, Node current, Node previous, vector<Node> subproblem_nodes, vector<vector<double>> *subproblem_requirements)
{
    auto current_id = graph.id(current);

    for (auto i=0; i < (int) subproblem_nodes.size(); i++)
    {
        auto node_id = graph.id(subproblem_nodes[i]);
        (*subproblem_requirements)[base_index][i] += requirements[node_id][current_id];
        (*subproblem_requirements)[i][base_index] += requirements[node_id][current_id];
    }

    for (auto e = Tree::IncEdgeIt(tree,current); e != INVALID; ++e)
    {
        auto node = tree.oppositeNode(current, e);
        if (node == previous)
        {
            continue;
        }
        addRequirements(base_index, node, current, subproblem_nodes, subproblem_requirements);
    }
}


vector<vector<double>> generateSubproblemsReq(const vector<Node>& subproblem_nodes)
{
    vector<vector<double>> subproblem_requirements = {};
    for (auto i=0; i < (int) subproblem_nodes.size(); i++)
    {
        subproblem_requirements.push_back(vector<double>());
        for (auto j=0; j < (int) subproblem_nodes.size(); j++)
        {
            if (j > i)
            {
                auto u_id = graph.id(subproblem_nodes[i]);
                auto v_id = graph.id(subproblem_nodes[j]);
                subproblem_requirements[i].push_back(requirements[u_id][v_id]);
            }
            else
            {
                subproblem_requirements[i].push_back(0);
            }
        }
    }

    for (auto i=0; i < (int) subproblem_nodes.size(); i++)
    {

        auto u = subproblem_nodes[i];
        // itera sobre as arestas,
        for (auto e = Tree::IncEdgeIt(tree,u); e != INVALID; ++e)
        {
            auto v = tree.oppositeNode(u, e);
            // se o vertice oposto não estiver no subproblema,
            // vai somando todos os requisitos dos outros pra esse
            if (find(subproblem_nodes.begin(), subproblem_nodes.end(), v) == subproblem_nodes.end())
            {
                addRequirements(i, v, u, subproblem_nodes, &subproblem_requirements);
            }
        }
    }

    for (auto i=0; i < (int) subproblem_nodes.size(); i++)
    {
        subproblem_requirements[i][i]=0;
    }

    return subproblem_requirements;
}


/*
    FORMULAÇÃO MATEMATICA COMEÇA AQUI
*/

void initVars(GRBModel& model, const vector<Node>& subproblem_nodes, 
              vector<uv>& pair_keys, vector<ouv>& triple_keys,
              map<uv, GRBVar>& x_map, map<ouv, GRBVar>& f_map, map<ouv, GRBVar>& y_map)
{
    // Creating vars
    for (int i=0; i < (int) subproblem_nodes.size(); i++)
    {
        for (int j=i; j < (int) subproblem_nodes.size(); j++)
        {
            auto u = subproblem_nodes[i];
            auto v = subproblem_nodes[j];
            auto e = findEdge(graph, u, v);
            if (e != INVALID)
            {
                // edges_tree[e] = false;
                int u_id = graph.id(u);
                int v_id = graph.id(v);
                // x_{o,u} (x for every edge (o,u) in subgraph)
                uv node_pair(u_id, v_id);
                pair_keys.push_back(node_pair);
                stringstream s;
                s << "x(" << u_id << "," << v_id << ")";
                auto x = model.addVar(0,1,0,GRB_BINARY, s.str());
                x_map.emplace(node_pair, x);
            }
        }
    }

    for (auto o : subproblem_nodes)
    {
        for (int i=0; i < (int) subproblem_nodes.size(); i++)
        {
            auto u = subproblem_nodes[i];
            for (int j=0; j < (int) subproblem_nodes.size(); j++)
            {
                auto v = subproblem_nodes[j];
                int o_id = graph.id(o);
                int u_id = graph.id(u);
                int v_id = graph.id(v);
                // By definition, we cannot have a (u,u) edge, so we don't need to test it
                if (findEdge(graph, u, v) != INVALID)
                {
                    // f^o_{u,v} for every arc in graph (u,v) != (v,u)
                    ouv node_triple(o_id, u_id, v_id);
                    triple_keys.push_back(node_triple);
                    stringstream s;
                    s << "f(" << o_id << "," << u_id << "," << v_id << ")";
                    auto f = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s.str());
                    f_map.emplace(node_triple, f);
                    stringstream s2;
                    s2 << "y(" << o_id << "," << u_id << "," << v_id << ")";
                    auto y = model.addVar(0,1,0,GRB_BINARY, s2.str());
                    y_map.emplace(node_triple, y);
                }
            }
        }
    }
}

void objective(GRBModel& model, const vector<Node>& subproblem_nodes, map<ouv, GRBVar>& f_map)
{
    GRBLinExpr objective_expr(0);
    for (int i=0; i < (int) subproblem_nodes.size(); i++)
    {
        auto u = subproblem_nodes[i];
        auto u_id = graph.id(u);
        for (int j=i; j < (int) subproblem_nodes.size(); j++)
        {
            auto v = subproblem_nodes[j];
            auto v_id = graph.id(v);
            auto e = findEdge(graph, u, v);
            if (e != INVALID)
            {
                // c_{u,v}
                auto length = lengths[e];
                GRBLinExpr inner_sum(0);
                for (auto o : subproblem_nodes)
                {
                    auto o_id = graph.id(o);
                    ouv triple_in(o_id,u_id,v_id);
                    ouv triple_out(o_id,v_id,u_id);
                    inner_sum += (f_map[triple_in] + f_map[triple_out]);
                }
                objective_expr += length * inner_sum;
            }
        }
    }
    model.setObjective(objective_expr, GRB_MINIMIZE);
}


// Origin sends the sum of all its requirements as initial flow
// OK
void first_constraint(GRBModel& model, const vector<Node>& subproblem_nodes, 
                      const vector<vector<double>>& subproblem_requirements,
                      map<ouv, GRBVar>& f_map)
{
    for (int i=0; i < (int) subproblem_nodes.size(); i++)
    {
        auto o = subproblem_nodes[i];
        int o_id = graph.id(o);

        GRBLinExpr left_sum(0);
        for (int j=0; j < (int) subproblem_nodes.size(); j++)
        {
            if (i==j) continue;

            auto u = subproblem_nodes[j];
            int u_id = graph.id(u);
            auto e = findEdge(graph, o, u);
            if (e != INVALID)
            {
                ouv triple(o_id,o_id,u_id);
                left_sum += f_map[triple];
            }
        }

        // O_o
        double right_sum=0;
        for (int j=0; j < (int) subproblem_nodes.size(); j++)
        {
            if (subproblem_requirements[i][j] >= 0)
            {
                right_sum += subproblem_requirements[i][j];
            }
        }
        model.addConstr(left_sum, GRB_EQUAL, right_sum);
    }
}


// Flow conservation
void second_constraint(GRBModel& model, const vector<Node>& subproblem_nodes, 
                       const vector<vector<double>>& subproblem_requirements, 
                       map<ouv, GRBVar>& f_map)
{
    // For every vertex
    for (int i=0; i < (int) subproblem_nodes.size(); i++)
    {
        auto o = subproblem_nodes[i];
        auto o_id = graph.id(o);
        // for every arc going from u
        for (int j=0; j < (int) subproblem_nodes.size(); j++)
        {
            if (i==j) 
            {
                continue;
            }

            auto u = subproblem_nodes[j];
            auto u_id = graph.id(u);

            GRBLinExpr sum_out(0); // Flow of o going out of u (f_{o,u,v})
            for (int k=0; k < (int) subproblem_nodes.size(); k++)
            {
                auto v = subproblem_nodes[k];
                // If the edge (u,v) exists, add f(o,u,v) to constraint
                if (findEdge(graph,u,v) != INVALID)
                {
                    int v_id = graph.id(v);
                    ouv triple_out(o_id,u_id,v_id);
                    sum_out += f_map[triple_out];
                }
            }

            GRBLinExpr sum_in(0); // Flow of o coming to u (f_{o,v,u})
            for (int k=0; k < (int) subproblem_nodes.size(); k++)
            {
                auto v = subproblem_nodes[k];
                auto v_id = graph.id(v);
                if (findEdge(graph,v,u) != INVALID)
                {
                    ouv triple_in(o_id,v_id,u_id);
                    sum_in += f_map[triple_in];
                }
            }
            model.addConstr(sum_out - sum_in, GRB_EQUAL, - subproblem_requirements[i][j]);
        }
    }
}

// f(o,u,v) must b zero if arc is not used
// Technically, it's the third and fourth constraint,
// but using a single loop in all nodes set (avoiding u==v), 
// it's possible to use the same constraint
void third_constraint(GRBModel& model, const vector<Node>& subproblem_nodes, 
                      vector<vector<double>>& subproblem_requirements, 
                      map<ouv, GRBVar>& f_map, map<ouv, GRBVar>& y_map)
{
    // For each origin
    for (int i=0; i < (int) subproblem_nodes.size(); i++)
    {
        auto o = subproblem_nodes[i];
        int o_id = graph.id(o);

        // That's the big-M
        // Sum of all requirements of o
        double big_m = 0;
        double min = INT32_MAX; // avoid NaN error using int max instead dbl max
        for (int j=0; j < (int) subproblem_nodes.size(); j++)
        {
            auto r = subproblem_requirements[i][j];
            if (r >= 0)
            {
                if (r < min)
                {
                    min = r;
                }
                big_m += r;
            }
        }
        if (min == INT32_MAX) min = 0;
        big_m -= min;

        // For each u
        for (int j=0; j < (int) subproblem_nodes.size(); j++)
        {
            auto u = subproblem_nodes[j];
            auto u_id = graph.id(u);
            // For each v
            // k=0 and only one constraint have the same result
            for (int k=i; k < (int) subproblem_nodes.size(); k++)
            {
                auto v = subproblem_nodes[k];
                auto e = findEdge(graph, u, v);
                if (e != INVALID)
                {
                    auto v_id = graph.id(v);
                    ouv triple_in(o_id,v_id,u_id);
                    ouv triple_out(o_id,u_id,v_id);
                    model.addConstr(f_map[triple_in], GRB_LESS_EQUAL, big_m*y_map[triple_in]);
                    model.addConstr(f_map[triple_out], GRB_LESS_EQUAL, big_m*y_map[triple_out]);
                }
            }
        }
    }
}

// Tree constraint for y
// Each origin o must have exact n-1 y_{o,u,v} activated
// OK
void fourth_constraint(GRBModel& model, const vector<Node>& subproblem_nodes, map<ouv, GRBVar>& y_map)
{
    for (auto o : subproblem_nodes)
    {
        auto o_id = graph.id(o);
        GRBLinExpr tree_constr_linexp(0);
        for (int i=0; i < (int) subproblem_nodes.size(); i++)
        {
            auto u = subproblem_nodes[i];
            auto u_id = graph.id(u);
            for (int j=0; j < (int) subproblem_nodes.size(); j++)
            {
                auto v = subproblem_nodes[j];   
                if (v == u || v == o) continue;
                auto v_id = graph.id(v);
                auto e = findEdge(graph, u, v);
                if (e != INVALID)
                {
                    ouv triple(o_id, u_id, v_id);
                    tree_constr_linexp += y_map[triple];
                }
            }
        }
        model.addConstr(tree_constr_linexp, GRB_EQUAL, (subproblem_nodes.size()-1));
    }
}

// Associate x_{u,v} to y_{o,u,v} xor y_{o,v,u}
// OK
void fifth_constraint(GRBModel& model, const vector<Node>& subproblem_nodes, vector<uv> pair_keys,
                      map<uv, GRBVar>& x_map ,map<ouv, GRBVar>& y_map)
{
    for (auto o : subproblem_nodes)
    {
        auto o_id = graph.id(o);
        // For each x_{i,j}: using pair_keys we mantain the v > u restriction
        for (auto pair : pair_keys)
        {
            int u_id = get<0>(pair);
            int v_id = get<1>(pair);
            ouv l(o_id,u_id,v_id);
            ouv r(o_id,v_id,u_id);
            model.addConstr(y_map[l] + y_map[r], GRB_LESS_EQUAL, x_map[pair]);
        }
    }
}


// Avoid cicle constraint
// Similar to fourth constraint
void sixth_constraint(GRBModel& model, const vector<Node>& subproblem_nodes,
                      vector<uv> pair_keys, map<uv, GRBVar>& x_map)
{
    GRBLinExpr cicle_constr_linexp(0);
    // constraint: sum(x) <= n-1
    for (auto pair : pair_keys)
    {
        cicle_constr_linexp += x_map[pair];
    }
    model.addConstr(cicle_constr_linexp, GRB_EQUAL, (subproblem_nodes.size()-1));
}

// TODO: Minimizar o número de laços, por enquanto é só pra garantir que funciona
void solveSubproblem(const vector<Node> subproblem_nodes)
{
    // Cria nova tabela de requisitos copiando os valores originais
    auto subproblem_requirements {generateSubproblemsReq(subproblem_nodes)};

    try
    {
        GRBModel model {GRBModel(env)};

        vector<uv> pair_keys {};
        map<uv, GRBVar> x_map;

        // f_{o,u,v} means: the flow of o going from u to v
        vector<ouv> triple_keys {};
        map<ouv, GRBVar> f_map;
        map<ouv, GRBVar> y_map;

        initVars(model, subproblem_nodes, pair_keys, triple_keys, x_map, f_map, y_map);

        // First of all, the objective
        objective(model, subproblem_nodes, f_map);

        // The origin sends the sum of all its requirements as initial flow
        first_constraint(model, subproblem_nodes, subproblem_requirements, f_map);

        // Flow conservation
        second_constraint(model, subproblem_nodes, subproblem_requirements, f_map);

        // f(o,u,v) must b zero if arc is not used
        third_constraint(model, subproblem_nodes, subproblem_requirements, f_map, y_map);

        // Tree constraint for y
        fourth_constraint(model, subproblem_nodes, y_map);

        // y to x
        fifth_constraint(model, subproblem_nodes, pair_keys, x_map, y_map);

        //Avoid cicle constraint
        sixth_constraint(model, subproblem_nodes, pair_keys, x_map);
        
        // model.write("wrong.lp");
        model.optimize();
        int status = model.get(GRB_IntAttr_Status);
        if( status == GRB_OPTIMAL)
        {
            // cout << "Subproblem result: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
            for (auto pair : pair_keys)
            {
                auto x = x_map[pair];
                auto x_name  = x.get(GRB_StringAttr_VarName);
                auto x_value = (bool) x.get(GRB_DoubleAttr_X);
                auto e = findEdge(graph, graph.nodeFromId(get<0>(pair)), graph.nodeFromId(get<1>(pair)));
                if (e != INVALID)
                {
                    edges_tree[e] = x_value;
                }
            }
            tree = Tree(graph, edges_tree);
            
            while (!clusters.empty())
            {
                clusters.pop_back();
            }
            for (auto n : nodes)
            {
                in_some_cluster[n] = false;
            }
            current_cluster = nullptr;
            // reinicia clusters
            // printEdgesTree();
            divideTree(root, INVALID);
        }
        else {
            cout << "fail" << endl;
        }

    }
    catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }
}


// TODO: Melhorar essa seleção (tá O(mn^4), da pra fazer melhor, mas assim to garantindo que não vai repetir)
vector<Node> selectTwoClusters()
{
    shuffle(clusters.begin(), clusters.end(), default_random_engine(seed++));
    vector<Node> final_cluster = {};
    for (auto node : clusters[0])
    {
        final_cluster.push_back(node);
    }
    for (int i=1; i < (int) clusters.size(); i++)
    {
        for (auto u : clusters[i])
        {
            for (auto v : clusters[0])
            {
                if (u == v || (findEdge(tree, u, v) != INVALID))
                {
                    for (auto node : clusters[i])
                    {
                        if (find(final_cluster.begin(), final_cluster.end(), node) == final_cluster.end())
                        {
                            final_cluster.push_back(node);
                        }
                    }
                    clusters.erase(clusters.begin()+i);
                    clusters.erase(clusters.begin());
                    clusters.push_back(final_cluster);
                    return final_cluster;
                }
            }
        }
    }
    return final_cluster;
}



int main(int argc, char* argv[])
{
    bool debug = false;
    if (argc != 4)
    {
        perror("usage: ./heuristic inputFile clusterSize iterNum\n");
        return 1;
    }

    srand(seed);
    env.set("LogFile", "heuristic_solver.log");
    env.set("OutputFlag", "0");
    env.start();

    // Inicialization
    // Read instance
    k = stoi(argv[2]);
    int iterNum = atoi(argv[3]);
    readInstance(argv[1]);

    root = nodes[rand()%n];
    cout << "instance read" << endl;
    generateInitialSolution();
    cout << "initial solution generated" << endl;
    tree = Tree(graph, edges_tree);
    divideTree(root, INVALID);
    cout << "clusterized" << endl;
    cout << "clusters created: " << clusters.size() << endl;
    // printEdgesTree();
    // printClusters();
    cout << calculateObjective() << endl;

    int five_percent = iterNum / 20;
    cout << "[";
    for (int i = 0; i < iterNum; i++)
    // while(clusters.size() > 1)
    {
        auto some_cluster = selectTwoClusters();
        if (debug)
        {
            for (auto node : some_cluster)
            {
                cout << graph.id(node) << " ";
            }
            cout << endl;
        }
        solveSubproblem(some_cluster);
        if ((i % five_percent) == 0)
        {
            cout << "|";
        }
    }
    cout << "]";
    // printTree(root, INVALID);
    cout << calculateObjective() << endl;
    // printTree(root, INVALID);
    // printEdgesTree();

	return 0;
}
