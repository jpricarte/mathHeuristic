#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <gurobi_c++.h>
#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>

#define K 3
#define C 2
#define MAX_ITER 10

using namespace std;
using namespace lemon;

// Defining types because I'm lazy
typedef ListDigraph Digraph;
typedef ListDigraph::Node Node;
typedef ListDigraph::Arc Arc;
typedef ListDigraph::ArcMap<double> ArcMap;
typedef ListDigraph::NodeMap<double> NodeMap;
typedef ListDigraph::NodeMap<int> IntNodeMap;
typedef ListDigraph::ArcMap<bool> BoolArcMap;
typedef FilterArcs<ListDigraph> SpanningTree;
typedef FilterNodes<SpanningTree, SpanningTree::NodeMap<bool>> Cluster;


// Graph elements
Digraph graph;
vector<Node> nodes = {};
vector<Arc> arcs = {};
ArcMap lengths(graph);
vector<vector<double>> requirements = {};
vector<vector<double>> distances = {};
int n, m;
// Tree elements
BoolArcMap arcInTree(graph, false);
SpanningTree tree(graph,arcInTree);
SpanningTree::ArcMap<double> treeLengths(tree);
Node treeRoot;
IntNodeMap nSubNodes(graph);

// Clusters elements
vector<vector<Node>> clusters = {};

// Solver Elements
GRBEnv env(true);

/**
 * @brief convert a instance file in multiple structs, all global variables
 * 
 * @param filename 
 */
void readInstance(char filename[])
{
    // Stream and temporary vars
    ifstream instanceFile(filename);
    string e1, e2, e3;
    uint16_t u, v;
    double value;

    // Getting graph size (n for nodes and m for arcs) and resizing everything
    instanceFile >> e1 >> e2;
    n = stoi(e1);
    m = 2*stoi(e2);
    graph.reserveNode(n);
    nodes.resize(n);
    requirements.resize(n, vector<double>(n));
    distances.resize(n, vector<double>(n, 0.0));
    graph.reserveArc(m);
    arcs.resize(m);

    // Adding nodes
    for (auto i=0; i<n; i++)
    {
        nodes[i] = graph.addNode();
    }
    // Adding arcs and length
    for (auto i=0; i<m; i+=2)
    {
        instanceFile >> e1 >> e2 >> e3;
        u = stoi(e1);
        v = stoi(e2);
        value = stod(e3);
        arcs[i] = graph.addArc(nodes[u], nodes[v]);
        arcs[i+1] = graph.addArc(nodes[v], nodes[u]);
        lengths[arcs[i]] = value;
        lengths[arcs[i+1]] = value;
        arcInTree[arcs[i]] = false;
        arcInTree[arcs[i+1]] = false;
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

/**
 * @brief Generate a initial solution, currently using dijkstra algorithm
 * @todo Create a lazyless solution
 * @return BoolEdgeMap* 
 */
void generateInitialSolution(Node root)
{ 
    Dijkstra<ListDigraph, ArcMap> dij(graph, lengths);
    dij.run(root);
    treeRoot = root;
    for(auto node : nodes)
    {
        if (node != root)
        {
            arcInTree[dij.predArc(node)] = true;
            treeLengths[dij.predArc(node)] = lengths[dij.predArc(node)];
        }
    }
}

int countSubNodes(Node root)
{

}

vector<Node> getNodesFromSubtree(Node root)
{

}

//TODO: Make it better
/*
    - Should return empty vector if size == K
    - Should avoid create cluster with size different than K
*/
vector<Node> createClusters(Node root = treeRoot, bool debug = false)
{
    
}


// Calculate objective (fun)
double calculateObjectiveFun()
{
    
}

void calculateSubProblem()
{
    
}


void printClusters()
{
    cout << "We have " << clusters.size() << " cluster" << endl;
    for (auto cluster : clusters)
    {
        for (auto node : cluster)
            cout << graph.id(node) << " - ";
        cout << endl;
    }
}

void printSubNodesAmount()
{
    for (auto node : nodes)
    {
        cout << graph.id(node) << ": " << nSubNodes[node] << endl;
    }
}

void printTree(Node root) {
    for (SpanningTree::OutArcIt a(tree, root); a!=INVALID; ++a)
    {
        cout << graph.id(root) << " " << graph.id(graph.target(a)) << endl;
        printTree(graph.target(a));
    }
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        perror("usage: ./heuristic inputFile root");
        return 1;
    }

    srand(0xc0ffee);

    env.set("LogFile", "solver.log");
    env.start();

    // Inicialization
    // Read instance
    readInstance(argv[1]);
    // Create initial solution
    int rootIndex = atoi(argv[2]);
    if (rootIndex >= n || rootIndex < 0) {
        perror("using as root a node that doesn't exists");
        return 1;
    }
    generateInitialSolution(nodes[rootIndex]);
    // cout << "initial solution created" << endl;
    //printTree(treeRoot);
    countSubNodes(treeRoot);
    //printSubNodesAmount();
    createClusters();
    // cout << "Clusters created" << endl;
    //printClusters();
    double objective = calculateObjectiveFun();
    cout << objective << endl;

    // Now, the funny part begins:
    // for (int _i; _i<MAX_ITER; _i++)
    // {
    //     // Choose two random clusters

    //     // use Gurobi to solve optimally

    //     // recalculate objective
    // }
	return 0;
}