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
typedef tuple<int, int> uv;
typedef tuple<int, int, int> ouv;


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

    for (auto i=0; i<subproblem_nodes.size(); i++)
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

// TODO: Minimizar o número de laços, por enquanto é só pra garantir que funciona
void solveSubproblem(vector<Node> subproblem_nodes)
{
    // Cria nova tabela de requisitos copiando os valores originais
    vector<vector<double>> subproblem_requirements = {};
    for (auto i=0; i<subproblem_nodes.size(); i++)
    {
        subproblem_requirements.push_back(vector<double>());
        for (auto j=0; j<subproblem_nodes.size(); j++)
        {
            auto u_id = graph.id(subproblem_nodes[i]);
            auto v_id = graph.id(subproblem_nodes[j]);
            subproblem_requirements[i].push_back(requirements[u_id][v_id]);
        }
    }
    // Após isso, soma os requisitos amarrados externamente ao vértice
    for (auto i=0; i<subproblem_nodes.size(); i++)
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

    for (auto i=0; i<subproblem_nodes.size(); i++)
    {
        subproblem_requirements[i][i]=0;
        for (auto j=0; j<subproblem_nodes.size(); j++)
        {
            cout << subproblem_requirements[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;

    /*
     * Solution based in Fernández, Luna-Mota, Hildenbrandt, Reinelt and Wiesberg
     * 2013 formulation
     * doi.org/10.1016/j.endm.2013.05.079
     */
    try
    {
        GRBModel model = GRBModel(env);

        vector<uv> pair_keys={};
        map<uv, GRBVar> x_map;

        vector<ouv> triple_keys={};
        map<ouv, GRBVar> f_map;
        map<ouv, GRBVar> y_map;

        // Creating vars
        for (int i=0; i<subproblem_nodes.size(); i++)
        {
            for (int j=i; j<subproblem_nodes.size(); j++)
            {
                Node u = subproblem_nodes[i];
                Node v = subproblem_nodes[j];
                Edge e = findEdge(graph, u, v);
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
                    //x.set(GRB_DoubleAttr_Start,edges_tree[e]);
                    x_map.emplace(node_pair, x);
                }
            }
        }

        for (auto o : subproblem_nodes)
        {
            for (auto u : subproblem_nodes)
            {
                for (auto v : subproblem_nodes)
                {
                    int o_id = graph.id(o);
                    int u_id = graph.id(u);
                    int v_id = graph.id(v);
                    if (findEdge(graph, u, v) != INVALID)
                    {
                        // f^o_{u,v} for every arc in graph (u,v) != (v,u)
                        ouv node_triple(o_id, u_id, v_id);
                        triple_keys.push_back(node_triple);
                        stringstream s;
                        s << "f(" << o_id << "," << u_id << "," << v_id << ")";
                        auto f = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER, s.str());
                        f_map.emplace(node_triple, f);

                        stringstream s2;
                        s2 << "y(" << o_id << "," << u_id << "," << v_id << ")";
                        auto y = model.addVar(0,1,0,GRB_BINARY, s2.str());
                        y_map.emplace(node_triple, y);
                    }
                }
            }
        }

        // The origin sends the sum of all its requirements as initial flow
        for (int i=0; i<subproblem_nodes.size(); i++)
        {
            auto o = subproblem_nodes[i];
            int o_id = graph.id(o);

            GRBLinExpr left_sum(0);
            for (auto w : subproblem_nodes)
            {
                int w_id = graph.id(w);
                auto e = findEdge(graph, o, w);
                if (e != INVALID)
                {
                    ouv triple(o_id,o_id,w_id);
                    left_sum += f_map[triple];
                }
            }

            double right_sum=0;
            for (int j=i; j<subproblem_nodes.size(); j++)
            {
                Node d = subproblem_nodes[j];
                int d_id = graph.id(d);
                if (subproblem_requirements[i][j] >= 0)
                {
                    right_sum += subproblem_requirements[i][j];
                    // cout << o_id  << " + " << d_id << " <= " << subproblem_requirements[i][j] << endl;
                }
            }
            model.addConstr(left_sum, GRB_EQUAL, right_sum);
        }

        // Flow conservation
        // for all o, for all v execpt o
        for (int i=0; i<subproblem_nodes.size(); i++)
        {
            auto o = subproblem_nodes[i];
            int o_id = graph.id(o);
            for (int j=0; j<subproblem_nodes.size(); j++)
            {
                Node u = subproblem_nodes[j];
                int u_id = graph.id(u);

                GRBLinExpr sum_in(0);
                GRBLinExpr sum_out(0);
                for (auto v : subproblem_nodes)
                {
                    if (v==o || v==u)
                    {
                        continue;
                    }
                    auto v_id = graph.id(v);
                    ouv triple_in(o_id, u_id, v_id);
                    sum_in += f_map[triple_in];
                    ouv triple_out(o_id, v_id, u_id);
                    sum_out += f_map[triple_out];
                }
                model.addConstr(sum_in - sum_out, GRB_EQUAL, subproblem_requirements[i][j]);
            }
        }

        for (auto o : subproblem_nodes)
        {
            auto o_id = graph.id(o);
            for (int i=0; i<subproblem_nodes.size(); i++)
            {
                auto u = subproblem_nodes[i];
                auto u_id = graph.id(u);
                for (int j=i; j<subproblem_nodes.size(); j++)
                {
                    auto req_sum = 0;
                    if (subproblem_requirements[i][j] >= 0)
                    {
                        req_sum += subproblem_requirements[i][j];
                    }
                    auto v = subproblem_nodes[j];
                    auto v_id = graph.id(v);
                    ouv triple_in (o_id, u_id, v_id);
                    ouv triple_out (o_id, u_id, v_id);
                    uv e(u_id, v_id);
                    model.addConstr(f_map[triple_in] - f_map[triple_out], GRB_LESS_EQUAL, req_sum*x_map(e));
                }
            }
        }

        // Finally, the objective
        GRBLinExpr objective_expr(0);
        for (auto triple : triple_keys)
        {
            Node u = graph.nodeFromId(get<1>(triple));
            Node v = graph.nodeFromId(get<2>(triple));
            Edge e = findEdge(graph, u, v);
            objective_expr += lengths[e] * f_map[triple];
        }

        model.setObjective(objective_expr, GRB_MINIMIZE);
        model.write("wrong.lp");
        model.optimize();
        int status = model.get(GRB_IntAttr_Status);
        if( status == GRB_OPTIMAL)
        {
            cout << "Success" << endl;
            for (auto pair : pair_keys)
            {
                auto x = x_map[pair];
                auto x_name  = x.get(GRB_StringAttr_VarName);
                auto x_value = (bool) x.get(GRB_DoubleAttr_X);
                auto e = findEdge(graph, graph.nodeFromId(get<0>(pair)),graph.nodeFromId(get<1>(pair)));
                if (e != INVALID)
                {
                    // if (edges_tree[e] !=  x_value)
                    // {
                    //     cout << graph.id(e) << " -> " << x_value << endl;
                    // }
                        // if (x_value)
                        //     cout << x_name << " - " << x_value << endl;
                    edges_tree[e] = x_value;
                }
            }
            tree = Tree(graph, edges_tree);
            // for (auto triple : triple_keys)
            // {
            //     auto f_name  = f_map[triple].get(GRB_StringAttr_VarName);
            //     auto f_value = f_map[triple].get(GRB_DoubleAttr_X);
            //     if (f_value != 0)
            //         cout << f_name << " - " << f_value << endl;
            // }
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
    random_shuffle(clusters.begin(), clusters.end());
    vector<Node> final_cluster = {};
    for (auto node : clusters[0])
    {
        final_cluster.push_back(node);
    }
    for (int i=1; i < clusters.size(); i++)
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
                    return final_cluster;
                }
            }
        }
    }
    return final_cluster;
}



int main(int argc, char* argv[])
{
    bool debug = true;
    if (argc != 4)
    {
        perror("usage: ./heuristic inputFile clusterSize iterNum\n");
        return 1;
    }

    srand(0xf0da5e);
    env.set("LogFile", "heuristic_solver.log");
    env.set("OutputFlag", "0");
    env.start();

    // Inicialization
    // Read instance
    k = stoi(argv[2]);
    readInstance(argv[1]);


    // for (auto linha : requirements)
    // {
    //     // for (auto elem : linha)
    //     // {
    //     //     cout << elem << "\t";
    //     // }
    //     cout << linha.size();
    //     cout << endl;
    // }

    // cout << "32, 2 - " << requirements[32][2] << endl;
    // cout << "6, 5 - " << requirements[6][5] << endl;

    root = nodes[rand()%n];
    cout << "instance read" << endl;
    generateInitialSolution();
    cout << "initial solution generated" << endl;
    tree = Tree(graph, edges_tree);
    divideTree(root, INVALID);
    cout << "clusterized" << endl;
    cout << "clusters created: " << clusters.size() << endl;
    // printClusters();
    cout << calculateObjective() << endl;
    for (int i = 0; i < atoi(argv[3]); i++)
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
            printEdgesTree();

    }
    // printTree(root, INVALID);
    cout << calculateObjective() << endl;
    // printTree(root, INVALID);
    printEdgesTree();

	return 0;
}
