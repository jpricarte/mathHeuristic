#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <random>
#include <gurobi_c++.h>
#include <lemon/list_graph.h>

#define K 3
#define C 2
#define MAX_ITER 10

using namespace std;
using namespace lemon;

// Defining types
typedef lemon::ListGraph Graph;
typedef lemon::ListGraph::Node Node;
typedef lemon::ListGraph::Edge Edge;

// Defining Graph elements
Graph graph;
int n, m; // Number of nodes and edges, respectively
unordered_set<Node> nodes;
unordered_set<Edge> edges;

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
    
	return 0;
}