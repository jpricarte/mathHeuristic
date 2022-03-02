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

    // Getting graph size (n for nodes and m for edges) and resizing everything
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
        for (auto j=i; j<n; j++)
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
        int counter = 0;
        while (!dij.emptyQueue())
        {
            counter++;
            Node v = dij.processNextNode();
            double distance = dij.dist(v);
            double requirement = requirements[tree.id(u)][tree.id(v)];
            cost += distance * requirement;
        }
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

// TODO: Minimizar o número de laços, por enquanto é só pra garantir que funciona
void solveSubproblem(vector<Node> subproblem_nodes)
{
    try 
    {
        GRBModel model = GRBModel(env);

        vector<uv> pair_keys={};
        map<uv, GRBVar> x_map;

        vector<ouv> triple_keys={};
        map<ouv, GRBVar> f_map;

        // Creating vars
        for (int i=0; i<subproblem_nodes.size(); i++)
        {
            for (int j=i; j<subproblem_nodes.size(); j++)
            {
                Node u = subproblem_nodes[i];
                Node v = subproblem_nodes[j];
                if (findEdge(graph, u, v) != INVALID)
                {
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
                    }
                }
            }
        }

        //Avoid cicle constraint
        GRBLinExpr cicle_constr_linexp(0);
        // First constraint: sum(x) <= n-1
        for (auto x : pair_keys)
        {
            cicle_constr_linexp += x_map[x];
        }
        model.addConstr(cicle_constr_linexp, GRB_EQUAL, (n-1));


        // Second constraint
        // origin doesn't receive it's own flow
        for (auto o : subproblem_nodes)
        {
            GRBLinExpr left_sum(0);
            int o_id = graph.id(o);
            for (auto u : subproblem_nodes)
            {
                int u_id = graph.id(u);
                auto e = findEdge(graph, u, o);
                if (e != INVALID)
                {
                    ouv triple(o_id,u_id,o_id);
                    left_sum += f_map[triple];
                }
            }

            model.addConstr(left_sum,GRB_EQUAL,0);
        }

        // Third constraint
        // Flow conservation
        // for all o, for all v execpt o
        for (auto o : subproblem_nodes)
        {
            int o_id = graph.id(o);
            for (auto v : subproblem_nodes)
            {
                // v cannot be equal to o
                if (v == o)
                {
                    continue;
                }
                int v_id = graph.id(v);

                GRBLinExpr sum_in(0);
                for (auto u : subproblem_nodes)
                {
                    if (findEdge(graph,u,v) != INVALID)
                    {
                        int u_id = graph.id(u);
                        ouv triple_in(o_id,u_id,v_id);
                        sum_in += f_map[triple_in];
                    }
                }

                GRBLinExpr sum_out(0);
                for (auto w : subproblem_nodes)
                {
                    if (findEdge(graph,v,w) != INVALID)
                    {
                        int w_id = graph.id(w);
                        ouv triple_out(o_id,v_id,w_id);
                        sum_out += f_map[triple_out];
                    }
                }
                model.addConstr(sum_in - sum_out, GRB_EQUAL, (int) requirements[o_id][v_id]);   
            }
        }
        
        // Fourth constraint
        // The origin sends the sum of all its requirements as initial flow
        for (auto o : subproblem_nodes)
        {
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
            for (auto d : subproblem_nodes)
            {
                int d_id = graph.id(d);
                if (o_id != d_id && requirements[o_id][d_id] > 0)
                {
                    right_sum += requirements[o_id][d_id];
                }
            }
            model.addConstr(left_sum, GRB_EQUAL, right_sum);
        }

        // Fifth constraint
        // flows from each origin can't be greater than the initial flow if the edge is used
        for (auto o : subproblem_nodes)
        {
            int o_id = graph.id(o);

            int right_sum = 0;
            for (auto d : subproblem_nodes)
            {
                int d_id = graph.id(d);
                if (o_id != d_id && requirements[o_id][d_id] > 0)
                {
                    right_sum += requirements[o_id][d_id];
                }
            }

            // For each edge (the pair ensures the edge exists)
            for (auto pair : pair_keys)
            {
                int u_id = get<0>(pair);
                int v_id = get<1>(pair);
                ouv l(o_id,u_id,v_id);
                ouv r(o_id,v_id,u_id);
                GRBLinExpr left_side(0);
                left_side += f_map[l];
                left_side += f_map[r];
                model.addConstr(left_side, GRB_LESS_EQUAL, right_sum * x_map[pair]);
            }
        }

        // Finally, the objective
        GRBLinExpr objective_expr(0);
        // o is the flow origin
        for (auto o : subproblem_nodes)
        {
            int o_id = graph.id(o);
            // left side of the arc
            for (auto u : subproblem_nodes)
            {
                int u_id = graph.id(u);
                // right side of the arc
                for (auto v : subproblem_nodes)
                {
                    int v_id = graph.id(v);
                    // Find if exists edge (u,v) (it means that exists the arcs (u,v) and (v,u))
                    Edge e = findEdge(graph, u, v);
                    if (e != INVALID)
                    {
                        ouv triple(o_id,u_id,v_id);
                        objective_expr += lengths[e] * f_map[triple];
                    }
                }
            }
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
                auto x_value = (bool) x.get(GRB_DoubleAttr_X);
                auto e = findEdge(graph, graph.nodeFromId(get<0>(pair)), graph.nodeFromId(get<1>(pair)));
                if (e != INVALID)
                {
                    edges_tree[e] = x_value;
                }
            }
            tree = Tree(graph, edges_tree);
        } else {
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


// TODO: Melhorar essa seleção
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
                        final_cluster.push_back(node);
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
        solveSubproblem(some_cluster);
        
    }
    // cout << calculateObjective() << endl;
    printTree(root, INVALID);
    
	return 0;
}