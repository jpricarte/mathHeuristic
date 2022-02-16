#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
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
typedef lemon::FilterEdges<Graph, vector<Node>> Subproblem;

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
GRBEnv env(true);

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
        while (!dij.emptyQueue())
        {
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

// void solveSubproblem(vector<Node> subproblem_nodes)
// {
//     Subproblem subproblem(graph, subproblem_nodes);
//     vector<Edge> subproblem_edges = {};
//     for (Subproblem::EdgeIt e(subproblem); e!=INVALID; ++e)
//     {
//         subproblem_edges.push_back(e);
//     }


//     GRBModel model(env);
//     // x_e means that e is part of the solution
//     int x_size = countEdges(subproblem);
//     GRBVar* x = model.addVars(x_size, GRB_BINARY);
//     // y_{uv}^e means that e is part of uv path
//     int y_size = (countNodes(subproblem)^2)*countEdges(subproblem);
//     GRBVar* y = model.addVars(y_size, GRB_BINARY);
//     // w_{uv}^w means that w is part of uv path
//     int w_size = (countNodes(subproblem))^3;
//     GRBVar* w = model.addVars(w_size, GRB_BINARY);
//     // Minimize the sum of all 
// }

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        perror("usage: ./heuristic inputFile root clusterSize\n");
        return 1;
    }

    srand(0xc0ffee);
    // env.set("LogFile", "solver.log");
    // env.start();

    // Inicialization
    // Read instance
    int root_index = stoi(argv[2]);
    k = stoi(argv[3]);
    readInstance(argv[1]);
    root = nodes[root_index];
    cout << "instance read" << endl;
    generateInitialSolution();
    cout << "initial solution generated" << endl;
    // printTree(root, INVALID);
    divideTree(root, INVALID);
    cout << "clusterized" << endl;
    cout << "clusters created: " << clusters.size();
    // printClusters();
    
	return 0;
}