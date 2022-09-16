//
// Created by jpricarte on 01/09/22.
//

#ifndef MATH_HEURISTIC_INSTANCE_H
#define MATH_HEURISTIC_INSTANCE_H
#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include <gurobi_c++.h>

typedef lemon::ListGraph Graph;
typedef lemon::ListGraph::Node Node;
typedef lemon::ListGraph::Edge Edge;

typedef lemon::ListGraph::EdgeMap<double> EdgeMapDouble;
typedef lemon::ListGraph::EdgeMap<bool> EdgeMapBool;

typedef lemon::FilterEdges<Graph,EdgeMapBool> Tree;
typedef lemon::FilterEdges<Graph,EdgeMapBool>::NodeMap<bool> TreeNodeMapBool;
typedef lemon::FilterNodes<Tree, TreeNodeMapBool> SubTree;

typedef std::tuple<int, int> uv;
typedef std::tuple<int, int, int> ouv;

struct Instance {
    unsigned int n;
    unsigned int m;
    Graph graph{};
    std::vector<Node> nodes{};
    std::vector<Edge> edges{};
    EdgeMapDouble lengths{graph};
    std::vector<std::vector<double>> requirements{};

    inline Instance() {
        n = 0;
        m = 0;
    };
    inline ~Instance() = default;
};

struct Solution{
public:
    Tree* tree;
    EdgeMapBool* edge_in_tree;
    double value;
    std::vector<std::vector<Node>> clusters;
    std::vector<Node>* current_cluster;

    inline Solution(Instance& i, EdgeMapBool& edge_in_tree, double value) :
                    value(value) {
        this->edge_in_tree = new EdgeMapBool{i.graph, false};
        for (auto edge : i.edges)
        {
            this->edge_in_tree->operator[](edge) = edge_in_tree[edge];
        }
        tree = new Tree{i.graph, *(this->edge_in_tree)};
        clusters = {};
        current_cluster = new std::vector<Node>{};
    };
public:
    void addToCluster(Node n, TreeNodeMapBool& in_some_cluster, const int cluster_size) {
        in_some_cluster[n] = true;
        current_cluster->push_back(n);

        // If the cluster is full, add cluster to clusters list
        if (current_cluster->size() == cluster_size)
        {
            clusters.push_back(*current_cluster);
            // printClusters();
            current_cluster = new std::vector<Node>{};
        }
    }
};

struct LinearProgram {
    GRBModel model;

    std::vector<uv> pair_keys;
    std::map<uv, GRBVar> x_map;

    // f_{o,u,v} means: the flow of o going from u to v
    std::vector<ouv> triple_keys;
    std::map<ouv, GRBVar> f_map;
    std::map<ouv, GRBVar> y_map;
};

#endif //MATH_HEURISTIC_INSTANCE_H
