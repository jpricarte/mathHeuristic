#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <map>
#include <tuple>
#include <random>
#include <algorithm>

#include <gurobi_c++.h>

#include <lemon/dfs.h>
#include <lemon/list_graph.h>
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
typedef lemon::FilterEdges<Graph,EdgeMapBool>::NodeMap<bool> TreeNodeMapBool;
typedef lemon::FilterNodes<Tree, TreeNodeMapBool> SubTree;
typedef tuple<int, int> uv;
typedef tuple<int, int, int> ouv;

const bool debug = false;
auto seed = 0xabeceda10;

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
EdgeMapBool edges_tree(graph, false);
Tree tree(graph, edges_tree);

// Cluster related things
TreeNodeMapBool in_some_cluster(tree, false);
vector<vector<Node>> clusters;
vector<Node>* current_cluster = nullptr;

// Solver Elements
GRBEnv env(false);


void readInstance(string filename)
{
    // Stream and temporary vars
    ifstream instanceFile(filename);
    string e1, e2, e3;
    unsigned int u, v;
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


// Print tree for debug
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

void printEdgesSubTree(SubTree t)
{
    for (SubTree::EdgeIt edge(t); edge != INVALID; ++edge)
    {
        if (edges_tree[edge])
        {
            cout << graph.id(graph.u(edge)) << " " << graph.id(graph.v(edge)) << endl;
        }
    }
}

// Print node clusters for debug
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

// Generator
double generateInitialSolutionDijkstra()
{
    for (auto edge : edges)
    {
        edges_tree[edge] = false;
    }
    Dijkstra<Graph, EdgeMapDouble> dij(graph, lengths);
    dij.run(root);
    double length_sum = 0.0;
    for (auto node : nodes)
    {
        Edge e = findEdge(graph, node, dij.predNode(node));
        if (e != INVALID)
        {
            if (!edges_tree[e]) length_sum += lengths[e];
            edges_tree[e] = true;
        }
    }
    return length_sum;
}

// Return true if the graph contains a cycle
bool containsCycle(Node n, Node up, TreeNodeMapBool& visited)
{
    bool has_cycle = false;
    visited[n] = true;
    for (Tree::IncEdgeIt e(tree, n); e != INVALID; ++e)
    {
        Node v = tree.oppositeNode(n,e);
        if (v != up)
        {
            if (visited[v] == true)
            {
                return true;
            }
            has_cycle = containsCycle(v, n, visited);
            // If some subgraph return a cycle, break the algorithm and return true
            if (has_cycle) return has_cycle;
        }
    }
    return false;
}

// Return true if the graph is a tree
bool verifyTree()
{
    TreeNodeMapBool visited(tree, false);
    bool contains_cycle = containsCycle(root, INVALID, visited);
    if (contains_cycle)
    {
        cout << "graph contains cycle" << endl;
        return !contains_cycle;
    }
    bool connected = true;
    for (auto node : nodes)
    {
        if (!visited[node])
        {
            connected = false;
            cout << "graph contains more than one connected component" << endl;
            break;
        }
    }

    // Is a tree if not contains cycle and is connected
    return (!contains_cycle) && connected;
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
        while (!dij.emptyQueue())
        {
            Node v = dij.processNextNode();
            // cout << graph.id(u) << " " << graph.id(v) << endl;
            double distance = dij.dist(v);
            double requirement = requirements[tree.id(u)][tree.id(v)];
            cost += distance * requirement;
        }
    }
    return cost/2;
}

double calculateSubproblemObjective(vector<Node> subproblem_nodes,
                                    vector<vector<double>> subproblem_requirements)
{
    TreeNodeMapBool subtreeNodes(tree, false);
    for (auto node : subproblem_nodes)
    {
        subtreeNodes[node] = true;
    }

    SubTree subtree(tree, subtreeNodes);
    // printEdgesSubTree(subtree);
    SubTree::EdgeMap<double> subproblem_lengths(subtree);
    for (auto u : subproblem_nodes)
    {
        for (auto v : subproblem_nodes)
        {
            auto e = findEdge(tree, u, v);
            if (e != INVALID)
            {
                subproblem_lengths[e] = lengths[e];
            }
        }
    }

    double cost = 0;
    for (auto u : subproblem_nodes)
    {
        Dijkstra<SubTree, SubTree::EdgeMap<double>> dij(subtree, subproblem_lengths);
        dij.init();
        dij.addSource(u);
        while (!dij.emptyQueue())
        {
            Node v = dij.processNextNode();
            if (u == v) continue;
            cout << graph.id(u) << " " << graph.id(v) << endl;
            double distance = dij.dist(v);
            int u_index = find(subproblem_nodes.begin(), subproblem_nodes.end(), u) 
                        - subproblem_nodes.begin();
            int v_index = find(subproblem_nodes.begin(), subproblem_nodes.end(), v)
                        - subproblem_nodes.begin();
            double requirement = subproblem_requirements[u_index][v_index];
            cost += distance * requirement;
        }
    }
    return cost;
}

// Aux function for tree division
void addToCluster(Node n)
{
    // If the current cluster is empty (doesn't exists yet) create it and add the node
    if (current_cluster == nullptr)
    {
        current_cluster = new vector<Node>();
    }

    in_some_cluster[n] = true;
    current_cluster->push_back(n);

    // If the cluster is full, add cluster to clusters list
    if (current_cluster->size() == k)
    {
        clusters.push_back(*current_cluster);
        // printClusters();
        current_cluster = nullptr;
    }
}

// Generate clusters (subtrees)
bool divideTree(Node n, Node up)
{
    if (n==INVALID) 
        return false;

    if (!in_some_cluster[n])
    {
        // cout << "Node " << tree.id(n) << " added for the first time" << endl;
        addToCluster(n);
    }

    for (Tree::EdgeIt it(tree); it != INVALID; ++it)
    {
        auto next_node = tree.oppositeNode(n, it);
        if (next_node == up) continue;
        
        // Recursive call to next node, return true if node exists in tree
        bool node_exists = divideTree(next_node, n);
        if (node_exists)
        {
            // If the current cluster still not full, and this node wasn't in this cluster, add to the current cluster
            if (current_cluster != nullptr &&
                find(current_cluster->begin(), current_cluster->end(), n) == current_cluster->end())
            {
                // cout << "Node " << tree.id(n) << " added for connection" << endl;
                addToCluster(n);
            }
        }
    }
    // If is the root node, add current cluster to the clusters list
    if (up == INVALID && current_cluster != nullptr)
    {
        clusters.push_back(*current_cluster);
        current_cluster = nullptr;
    }
    return true;
}

Node findCentroidOfTree(Node n, SubTree& subtree, vector<Node> subproblem_nodes)
{
    Node min_node = n;
    int min_dist = INT32_MAX;
    for (auto u : subproblem_nodes)
    {
        int distance = 0;
        Dfs<SubTree> dfs(subtree);
        dfs.run(u);
        for (auto v : subproblem_nodes)
        {
            distance += dfs.dist(v);
        }
        if (distance <= min_dist)
        {
            min_dist = distance;
            min_node = u;
        }
    }
    return min_node;
}

void countNodesFromEdge(Node u, Node from, SubTree& subtree,
                SubTree::EdgeMap<int>& nodes_from_edge)
{
    Edge curr_edge = findEdge(subtree, u, from);
    for (SubTree::IncEdgeIt it(subtree,u); it != INVALID; ++it)
    {
        Node v = subtree.oppositeNode(u, it);
        if (v == from) continue;
        countNodesFromEdge(v, u, subtree, nodes_from_edge);
        if (curr_edge == INVALID) continue;
        nodes_from_edge[curr_edge] += nodes_from_edge[it];
    }
}

void clusterFromTreeByEdge(Node n, Edge last_edge, SubTree& subtree, vector<Node>& cluster)
{
    cluster.push_back(n);
    for (SubTree::IncEdgeIt it(subtree,n); it != INVALID; ++it)
    {
        if (it == last_edge) continue;
        Node u = subtree.oppositeNode(n,it);
        clusterFromTreeByEdge(u, it, subtree, cluster);
    }
}

// Split tree in two if contains a perfect centroid (tree_size/2)
bool splitInEdge(int tree_size, SubTree& subtree, SubTree::EdgeMap<int>& nodes_from_edge, int i, int j)
{
    // Iterate over edges to find a edge centroid
    for (SubTree::EdgeIt it(subtree); it != INVALID; ++it)
    {
        if (nodes_from_edge[it] == ceil(tree_size / 2))
        {
            // Divide in this edge
            Node u = subtree.u(it);
            Node v = subtree.v(it);
            vector<Node> cluster_u{};
            vector<Node> cluster_v{};
            clusterFromTreeByEdge(u, it, subtree, cluster_u);
            clusterFromTreeByEdge(v, it, subtree, cluster_v);
            clusters[i] = cluster_u;
            clusters[j] = cluster_v;
            return true;
        }
    }
    return false;
}

bool splitInNode(Node centroid, SubTree& subtree, int i, int j)
{
    // If contains more than 2 children, cannot use this approach
    int count=0;
    for (SubTree::IncEdgeIt it(subtree,centroid); it != INVALID; ++it)
        count++;
    if (count != 2) return false;
    // Iterate over edges to find a edge centroid
    vector<vector<Node>> new_clusters = {};
    for (SubTree::IncEdgeIt it(subtree, centroid); it != INVALID; ++it)
    {
            // Divide in this edge
            Node u = subtree.oppositeNode(centroid, it);
            vector<Node> cluster{centroid};
            clusterFromTreeByEdge(u, it, subtree, cluster);
            new_clusters.push_back(cluster);
    }
    clusters[i] = new_clusters[i];
    clusters[j] = new_clusters[j];
    return true;
}


void knapsackApproach(Node centroid, SubTree& subtree, SubTree::EdgeMap<int>& nodes_from_edge, int tree_size, int i, int j)
{
    GRBModel model_kp {GRBModel(env)};
    vector<GRBVar> vars;
    GRBLinExpr objective_expr(0);
    GRBLinExpr constraint_expr(0);

    for (SubTree::IncEdgeIt it(subtree,centroid); it != INVALID; ++it)
    {
        auto x = model_kp.addVar(0,1,0,GRB_BINARY);
        vars.push_back(x);
        objective_expr += x;
        constraint_expr += (nodes_from_edge[it] * x);
    }

    model_kp.addConstr(constraint_expr, GRB_LESS_EQUAL, (ceil(tree_size / 2)+1));
    model_kp.setObjective(objective_expr, GRB_MAXIMIZE);
    model_kp.optimize();
    int status = model_kp.get(GRB_IntAttr_Status);
    if( status == GRB_OPTIMAL)
    {
        int k=0;
        vector<Node> cluster_one{centroid};
        vector<Node> cluster_two{centroid};
        for (SubTree::IncEdgeIt it(subtree,centroid); it != INVALID; ++it)
        {
            Node u = subtree.oppositeNode(centroid, it);
            if (u == INVALID) continue;

            if (vars[k].get(GRB_DoubleAttr_X))
            {
                // Add all nodes bellow edge in a cluster
                clusterFromTreeByEdge(u, it, subtree, cluster_one);
            }
            else 
            {
                clusterFromTreeByEdge(u, it, subtree, cluster_two);
            }
            k++;
        }
        clusters[i] = cluster_one;
        clusters[j] = cluster_two;
        // cout << "new1 = " << cluster_one.size() << "; new2 = " << cluster_two.size() << endl;
        // if (abs(((int)cluster_one.size()) - ((int)cluster_two.size())) > 2 )
        // {
        //     cout << "centroid: " << graph.id(centroid) << "; cluster lengths: " << ceil(tree_size/2) << endl;
        //     printEdgesSubTree(subtree);
        // }
        // cout << "==========" << endl;
    }
}

void splitTree(SubTree& subtree, vector<Node>& subproblem_nodes, int i, int j)
{
    // first of all, find a centroid
    Node centroid = findCentroidOfTree(subproblem_nodes[0], subtree, subproblem_nodes);
    
    // After, count nodes hanged in each edge
    SubTree::EdgeMap<int> nodes_from_edge(subtree, 1);
    countNodesFromEdge(centroid, INVALID, subtree, nodes_from_edge);
    
    // After that, try to split the tree in an edge, cause it's simple
    bool divided = splitInEdge(subproblem_nodes.size(), subtree, nodes_from_edge, i, j);
    if (divided)
    {
        if (debug)
            cout << "divided in edge" << endl;
        return;
    }

    // Also, try to split in Node
    // divided = splitInNode(centroid, subtree, i, j);
    if (divided)
    {
        if (debug)
            cout << "divided in centroid" << endl;
        return;
    }

    // If nothing works, use a knapsack approach
    knapsackApproach(centroid, subtree, nodes_from_edge, subproblem_nodes.size(), i, j);
    if (debug)
        cout << "divided using knapsack approach" << endl;

}

// Aux function for subproblem requirements
/*
    base-index: indice no vetor do subproblema da aresta
a árvore que será adicionada
    current: vértice que está sendo analizado nesse momento
    previous: último vértice analisado, apenas para evitar volta
    outros: vetores comuns
*/
void addRequirements(int base_index, Node current , Node previous,
                     vector<Node> subproblem_nodes,
                     vector<vector<double>> *subproblem_requirements)
{
    auto current_id = graph.id(current);

    // Para cada vértice da árvore
    for (auto i=0; i < (int) subproblem_nodes.size(); i++)
    {
        // Inicialmente, o requerimento entre o vértice da árvore que será sobrecarregado (base_index)
        // e os outros vértices receberá uma soma do requerimento do novo vértice (current)
        auto node_id = graph.id(subproblem_nodes[i]);
        if (i > base_index)
        {
            (*subproblem_requirements)[base_index][i] += requirements[node_id][current_id];
        }
        else
        {
            (*subproblem_requirements)[i][base_index] += requirements[node_id][current_id];
        }
    }

    // Recursivamente, vai para os outros vértices ligados a esse que ainda não foram acessados
    // Por ser uma árvore, não precisamos de controle de acessados
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

        // Cria matriz auxiliar de requerimentos (triangular superior)
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
        // Para cada vértice u
        auto u = subproblem_nodes[i];
        // itera sobre as arestas,
        for (auto e = Tree::IncEdgeIt(tree,u); e != INVALID; ++e)
        {
            auto v = tree.oppositeNode(u, e);
            // se o vertice oposto não estiver no subproblema,
            // vai somando todos os requisitos dos vértices agregados para ele
            if (find(subproblem_nodes.begin(), subproblem_nodes.end(), v) == subproblem_nodes.end())
            {
                addRequirements(i, v, u, subproblem_nodes, &subproblem_requirements);
            }
        }
    }

    for (auto i=0; i < (int) subproblem_nodes.size(); i++)
    {
        subproblem_requirements[i][i] = 0;
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
                // x.set(GRB_DoubleAttr_Start, 1.0);
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

double solveSubproblem(vector<Node> subproblem_nodes, bool* was_modified, int i, int j)
{
    // Cria nova tabela de requisitos copiando os valores originais
    auto subproblem_requirements {generateSubproblemsReq(subproblem_nodes)};

    double init_value = calculateSubproblemObjective(subproblem_nodes, subproblem_requirements);

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
            auto solver_result = model.get(GRB_DoubleAttr_ObjVal);

            if (!verifyTree())
            {
                if (debug)
                    perror("new solution not valid!\n");
                return 0;
            }
            if (((int)solver_result) < ((int) init_value)) // Update tree only if we got a better solution
            {
                if (debug)
                {
                    cout << "optimized!" << endl;
                    cout << "Subproblem result: " << solver_result;
                    cout << "; Original result: " << init_value << endl;
                }
                for (auto pair : pair_keys)
                {
                    auto x = x_map[pair];
                    auto x_name  = x.get(GRB_StringAttr_VarName);
                    auto x_value = (bool) x.get(GRB_DoubleAttr_X);
                    auto u = graph.nodeFromId(get<0>(pair));
                    auto v = graph.nodeFromId(get<1>(pair));
                    auto e = findEdge(graph, u, v);
                    if (e != INVALID)
                    {
                        if (edges_tree[e] != (bool) x_value)
                        {
                            *was_modified = true;
                        }
                        edges_tree[e] = x_value;
                    }
                }

                tree = Tree(graph, edges_tree);
                TreeNodeMapBool in_subtree(tree, false);

                for (auto n : subproblem_nodes)
                {
                    in_subtree[n] = true;
                }
                SubTree subtree(tree, in_subtree);
                auto calcted_value = 1;// calculateSubproblemObjective(subproblem_nodes, subproblem_requirements);
                // printEdgesSubTree(subtree);
                splitTree(subtree, subproblem_nodes, i, j);

                // if (solver_result != calcted_value)
                // {
                    cout << "initial: " << init_value << endl;
                    cout << "c1 = " << i << "; c2 = " << j << endl;
                    cout << "solver: " << solver_result << "\tcalculated: " << calcted_value << endl;
                    printClusters();
                    // printEdgesTree();
                    cout << "----------------------------------------------------------" << endl;
                // }

                return solver_result - init_value;
            }
            else
            {
                if (debug)
                {
                    cout << "no optimized" << endl
                    << "Subproblem result: " << solver_result 
                    << "; Original result: " << init_value << endl;
                }
                return 0;
            }
        }
        else {
            if (debug)
                cout << "fail" << endl;
            return 0;
        }

    }
    catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    return 0;
}


// First way to select clusters
vector<Node> selectTwoClusters(int sel)
{
    // selected clusters
    vector<Node> final_cluster = {};
    // deep copy of first cluster
    for (auto node : clusters[sel])
    {
        final_cluster.push_back(node);
    }

    int i = rand() % final_cluster.size();
    Node sel_node = final_cluster[i];

    for (auto cluster : clusters)
    {
        if (cluster == clusters[sel])
            continue;

        for (auto node : cluster)
        {
            if (node == sel_node || findEdge(tree, node, sel_node) != INVALID)
            {
                for (auto node : cluster)
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

    return final_cluster;
}

// Second way to select clusters
vector<Node> selectTwoClusters(int first, int second)
{
    vector<Node> final_cluster = {};
    
    for (auto node_one : clusters[first])
    {
        for (auto node_two : clusters[second])
        {
            // If clusters are connected, create a merged cluster
            if (node_one == node_two || findEdge(tree, node_one, node_two) != INVALID)
            {
                for (auto node : clusters[first])
                {
                    final_cluster.push_back(node);
                }

                for (auto node : clusters[second])
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
    return final_cluster;
}

// Select root using min_paths approach
void chooseRoot()
{
    auto min_root = nodes[0];
    double min_length = INT32_MAX;
    for (auto node : nodes)
    {
        root = node;
        auto curr_length = generateInitialSolutionDijkstra();
        if (curr_length < min_length)
        {
            min_root = node;
            min_length = curr_length;
        }
    }
    root = min_root;
    generateInitialSolutionDijkstra();
}

int main(int argc, char* argv[])
{
    cout << "received ";
        for (int i=0; i<argc; i++)
        {
            cout << argv[i] << " ";
        }
        cout << endl;
    if (argc != 4)
    {
        perror("usage: ./heuristic inputFile clusterSize maxIter mode \n");
        return 1;
    }

    srand(seed);
    // env.set("LogFile", "heuristic_solver.log");
    env.set("OutputFlag", "0");
    env.start();

    // Inicialization
    // Read instance
    k = stoi(argv[2]);
    int max_iter = atoi(argv[3]);
    string filename(argv[1]);
    readInstance(filename);
    chooseRoot();
    cout << "instance read" << endl;
    generateInitialSolutionDijkstra();
    if (!verifyTree())
    {
        perror("Initial solution not valid!\n");
        printEdgesTree();
        return 1;
    }
    cout << "initial solution generated" << endl;
    tree = Tree(graph, edges_tree);
    divideTree(root, INVALID);
    cout << "clusterized" << endl;
    cout << "clusters created: " << clusters.size() << endl;
    double objective = calculateObjective();
    cout << objective << endl;
    auto start = chrono::high_resolution_clock::now();

    int iterNum = 0;
    vector<int> modified_clusters = {};

    for (int i = 0; i < clusters.size(); i++)
    {
        for (int j = i+1; j < clusters.size(); j++)
        {
            auto some_cluster = selectTwoClusters(i, j);
            if (some_cluster.empty())
                continue;
            iterNum++;

            bool was_modified = false;

            // return the diference between the original value and the new solution
            auto diff = solveSubproblem(some_cluster, &was_modified, i, j); 

            auto new_objective = objective + diff;
            if (new_objective < objective)
            {
                if (debug)
                    cout << "melhorou" << endl;
                objective = new_objective;
            }                

            // Add to modified_clusters if not in there yet
            if (was_modified)
            {
                if (find(modified_clusters.begin(), modified_clusters.end(), i) == modified_clusters.end())
                {
                    modified_clusters.push_back(i);
                }
                if (find(modified_clusters.begin(), modified_clusters.end(), j) == modified_clusters.end())
                {
                    modified_clusters.push_back(j);
                }
            }

            if (debug)
            {
                for (auto node : some_cluster)
                {
                    cout << graph.id(node) << " ";
                }
                cout << endl;
            }
        } 
    }

    // After run in every cluster at least once,
    // iterate until clusters dosen't modify anymore or until
    // the iteration limit is reached
    while (!modified_clusters.empty() && iterNum < max_iter)
    {
        int i = modified_clusters.front();
        modified_clusters.erase(modified_clusters.begin());
        for (int j = 0; j < clusters.size(); j++)
        {
            if (i == j) continue;

            if (debug)
            {
                cout << modified_clusters.size() << endl;
            }
            
            auto some_cluster = selectTwoClusters(i, j);
            if (some_cluster.empty())
                continue;
            iterNum++;

            bool was_modified = false;

            // return the diference between the original value and the new solution
            auto diff = solveSubproblem(some_cluster, &was_modified, i, j); 

            auto new_objective = objective + diff;
            if (new_objective < objective)
            {
                if (debug)
                    cout << "melhorou" << endl;
                objective = new_objective;
            }

            if (was_modified)
            {
                if (find(modified_clusters.begin(), modified_clusters.end(), i) == modified_clusters.end())
                {
                    modified_clusters.push_back(i);
                }
                if (find(modified_clusters.begin(), modified_clusters.end(), j) == modified_clusters.end())
                {
                    modified_clusters.push_back(j);
                }
            }


            if (debug)
            {
                for (auto node : some_cluster)
                {
                    cout << graph.id(node) << " ";
                }
                cout << endl;
            }
        }

    }

    auto stop = chrono::high_resolution_clock::now();
    auto exec_time = chrono::duration_cast<chrono::seconds>(stop - start);
    auto value = calculateObjective();
    printClusters();
    printEdgesTree();
    cout << "calculado: " << value << endl;
    cout << "atualizado: " << objective << endl;
    auto name_index = filename.find_last_of("/");
    auto ext_index = filename.find_last_of(".");
    auto stats_file = "./output/" + 
                    filename.substr(name_index, ext_index) + 
                    "_stats.csv";
    std::ofstream outfile;
    outfile.open(stats_file, std::ios_base::app); // append instead of overwrite
    outfile << k << "," << iterNum << "," << exec_time.count() << "," << value << endl;
    outfile.close();

	return 0;
}
