//
// Created by jpricarte on 01/09/22.
//

#ifndef MATH_HEURISTIC_HEURISTIC_H
#define MATH_HEURISTIC_HEURISTIC_H

#include "instance.h"
#include "auxiliar.h"
#include <lemon/dijkstra.h>
#include <gurobi_c++.h>
#include <lemon/dfs.h>
#include <cmath>


/**
 @name: generateIntialSolutionDijkstra
 @param: instance a Instance struct.
 @return: Solution - A Solution struct with the generated tree and the sum of lengths
 @brief: uses Dijkstra to generate a tree.
 @warning: It not returns the communication cost!
 * */
Solution* generateInitialSolutionDijkstra(Instance& instance, const Node root, int cluster_size)
{
    EdgeMapBool edges_tree(instance.graph, false);

    lemon::Dijkstra<Graph, EdgeMapDouble> dij(instance.graph, instance.lengths);
    dij.run(root);
    double length_sum = 0.0;
    for (auto node : instance.nodes)
    {

        Edge e = findEdge(instance.graph, node, dij.predNode(node));
        if (e != lemon::INVALID)
        {
            if (!edges_tree[e]) length_sum += instance.lengths[e];
            edges_tree[e] = true;
        }
    }

    return new Solution {instance, edges_tree, length_sum};
}

/**
 @name: chooseRoot
 @param: instance a Instance struct.
 @return: Node root, the selected root
 @brief: choose a root that gives the MST of an instance.
 * */
Node chooseRoot(Instance& instance, int cluster_size)
{
    Node min_root = instance.nodes[0];
    double min_length = INT32_MAX;
    for (auto node : instance.nodes)
    {
        auto curr_solution = generateInitialSolutionDijkstra(instance, node, cluster_size);
        if (curr_solution->value < min_length)
        {
            min_root = node;
            min_length = curr_solution->value;
        }
    }
    return min_root;
}

bool divideTree(Node n, Node up, TreeNodeMapBool& in_some_cluster, Solution* solution, int cluster_size)
{
    if (n==lemon::INVALID)
        return false;

    if (!in_some_cluster[n])
    {
        // cout << "Node " << tree.id(n) << " added for the first time" << endl;
        solution->addToCluster(n, in_some_cluster, cluster_size);
    }
    for (Tree::EdgeIt it(*solution->tree); it != lemon::INVALID; ++it)
    {
        auto next_node = solution->tree->oppositeNode(n, it);
        if (next_node == up) continue;

        // Recursive call to next node, return true if node exists in tree
        bool node_exists = divideTree(next_node, n, in_some_cluster, solution, cluster_size);
        if (node_exists)
        {
            // If the current cluster still not full, and this node wasn't in this cluster, add to the current cluster
            if (!solution->current_cluster->empty() &&
                find(solution->current_cluster->begin(), solution->current_cluster->end(), n)
                     == solution->current_cluster->end())
            {
                // cout << "Node " << tree.id(n) << " added for connection" << endl;
                solution->addToCluster(n, in_some_cluster, cluster_size);
            }
        }
    }
    // If is the root node, add current cluster to the clusters list
    if (up == lemon::INVALID && solution->current_cluster != nullptr)
    {
        solution->clusters.push_back(*solution->current_cluster);
        solution->current_cluster = new std::vector<Node>{};
    }

    return true;
}

double calculateObjective(Instance& instance, Solution* solution)
{
    double cost = 0;
    for (auto u : instance.nodes)
    {
        lemon::Dijkstra<Tree, EdgeMapDouble> dij(*solution->tree, instance.lengths);
        dij.init();
        dij.addSource(u);
        // Avoiding processing the u-u edge
        dij.processNextNode();
        while (!dij.emptyQueue())
        {
            Node v = dij.processNextNode();
            // cout << graph.id(u) << " " << graph.id(v) << endl;
            double distance = dij.dist(v);
            double requirement = instance.requirements[Graph::id(u)][Graph::id(v)];
            cost += distance * requirement;
        }
    }
    return cost/2;
}

double calculateSubproblemObjective(std::vector<Node> subproblem_nodes,
                                    std::vector<std::vector<double>> subproblem_requirements,
                                    Instance& instance, Solution* solution)
{
    TreeNodeMapBool subtreeNodes(*solution->tree, false);
    for (auto node : subproblem_nodes)
    {
        subtreeNodes[node] = true;
    }

    SubTree subtree(*solution->tree, subtreeNodes);
    // printEdgesSubTree(subtree);
    SubTree::EdgeMap<double> subproblem_lengths(subtree);
    for (auto u : subproblem_nodes)
    {
        for (auto v : subproblem_nodes)
        {
            auto e = findEdge(*solution->tree, u, v);
            if (e != lemon::INVALID)
            {
                subproblem_lengths[e] = instance.lengths[e];
            }
        }
    }

    double cost = 0;
    for (auto u : subproblem_nodes)
    {
        lemon::Dijkstra<SubTree, SubTree::EdgeMap<double>> dij(subtree, subproblem_lengths);
        dij.init();
        dij.addSource(u);
        while (!dij.emptyQueue())
        {
            Node v = dij.processNextNode();
            if (u == v) continue;
//            std::cout << Graph::id(u) << " " << Graph::id(v) << std::endl;
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

void initVars(GRBModel& model, const std::vector<Node>& subproblem_nodes,
              std::vector<uv>& pair_keys, std::vector<ouv>& triple_keys,
              std::map<uv, GRBVar>& x_map, std::map<ouv, GRBVar>& f_map, std::map<ouv, GRBVar>& y_map,
              const Instance& instance) {
    // Creating vars
    for (int i = 0; i < (int) subproblem_nodes.size(); i++) {
        for (int j = i; j < (int) subproblem_nodes.size(); j++) {
            auto u = subproblem_nodes[i];
            auto v = subproblem_nodes[j];
            auto e = findEdge(instance.graph, u, v);
            if (e != lemon::INVALID) {
                // edges_tree[e] = false;
                int u_id = Graph::id(u);
                int v_id = Graph::id(v);
                // x_{o,u} (x for every edge (o,u) in subgraph)
                uv node_pair(u_id, v_id);
                pair_keys.push_back(node_pair);
                std::stringstream s;
                s << "x(" << u_id << "," << v_id << ")";
                auto x = model.addVar(0, 1, 0, GRB_BINARY, s.str());
                // x.set(GRB_DoubleAttr_Start, 1.0);
                x_map.emplace(node_pair, x);
            }
        }
    }
}

void objective(GRBModel& model, const std::vector<Node>& subproblem_nodes, std::map<ouv, GRBVar>& f_map,
               const Instance& instance)
{
    GRBLinExpr objective_expr(0);
    for (int i=0; i < (int) subproblem_nodes.size(); i++)
    {
        auto u = subproblem_nodes[i];
        auto u_id = Graph::id(u);
        for (int j=i; j < (int) subproblem_nodes.size(); j++)
        {
            auto v = subproblem_nodes[j];
            auto v_id = Graph::id(v);
            auto e = findEdge(instance.graph, u, v);
            if (e != lemon::INVALID)
            {
                // c_{u,v}
                auto length = instance.lengths[e];
                GRBLinExpr inner_sum(0);
                for (auto o : subproblem_nodes)
                {
                    auto o_id = Graph::id(o);
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
void first_constraint(GRBModel& model, const std::vector<Node>& subproblem_nodes,
                      const std::vector<std::vector<double>>& subproblem_requirements,
                      std::map<ouv, GRBVar>& f_map, const Instance& instance)
{
    for (int i=0; i < (int) subproblem_nodes.size(); i++)
    {
        auto o = subproblem_nodes[i];
        int o_id = Graph::id(o);

        GRBLinExpr left_sum(0);
        for (int j=0; j < (int) subproblem_nodes.size(); j++)
        {
            if (i==j) continue;

            auto u = subproblem_nodes[j];
            int u_id = Graph::id(u);
            auto e = findEdge(instance.graph, o, u);
            if (e != lemon::INVALID)
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
void second_constraint(GRBModel& model, const std::vector<Node>& subproblem_nodes,
                       const std::vector<std::vector<double>>& subproblem_requirements,
                       std::map<ouv, GRBVar>& f_map, const Instance& instance)
{
    // For every vertex
    for (int i=0; i < (int) subproblem_nodes.size(); i++)
    {
        auto o = subproblem_nodes[i];
        auto o_id = Graph::id(o);
        // for every arc going from u
        for (int j=0; j < (int) subproblem_nodes.size(); j++)
        {
            if (i==j)
            {
                continue;
            }

            auto u = subproblem_nodes[j];
            auto u_id = Graph::id(u);

            GRBLinExpr sum_out(0); // Flow of o going out of u (f_{o,u,v})
            for (int k=0; k < (int) subproblem_nodes.size(); k++)
            {
                auto v = subproblem_nodes[k];
                // If the edge (u,v) exists, add f(o,u,v) to constraint
                if (findEdge(instance.graph,u,v) != lemon::INVALID)
                {
                    int v_id = Graph::id(v);
                    ouv triple_out(o_id,u_id,v_id);
                    sum_out += f_map[triple_out];
                }
            }

            GRBLinExpr sum_in(0); // Flow of o coming to u (f_{o,v,u})
            for (int k=0; k < (int) subproblem_nodes.size(); k++)
            {
                auto v = subproblem_nodes[k];
                auto v_id = Graph::id(v);
                if (findEdge(instance.graph,v,u) != lemon::INVALID)
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
void third_constraint(GRBModel& model, const std::vector<Node>& subproblem_nodes,
                      std::vector<std::vector<double>>& subproblem_requirements,
                      std::map<ouv, GRBVar>& f_map, std::map<ouv, GRBVar>& y_map,
                      const Instance& instance)
{
    // For each origin
    for (int i=0; i < (int) subproblem_nodes.size(); i++)
    {
        auto o = subproblem_nodes[i];
        int o_id = Graph::id(o);

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
            auto u_id = Graph::id(u);
            // For each v
            // k=0 and only one constraint have the same result
            for (int k=i; k < (int) subproblem_nodes.size(); k++)
            {
                auto v = subproblem_nodes[k];
                auto e = findEdge(instance.graph, u, v);
                if (e != lemon::INVALID)
                {
                    auto v_id = Graph::id(v);
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
void fourth_constraint(GRBModel& model, const std::vector<Node>& subproblem_nodes,
                       std::map<ouv, GRBVar>& y_map, const Instance& instance)
{
    for (auto o : subproblem_nodes)
    {
        auto o_id = Graph::id(o);
        GRBLinExpr tree_constr_linexp(0);
        for (int i=0; i < (int) subproblem_nodes.size(); i++)
        {
            auto u = subproblem_nodes[i];
            auto u_id = Graph::id(u);
            for (int j=0; j < (int) subproblem_nodes.size(); j++)
            {
                auto v = subproblem_nodes[j];
                if (v == u || v == o) continue;
                auto v_id = Graph::id(v);
                auto e = findEdge(instance.graph, u, v);
                if (e != lemon::INVALID)
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
void fifth_constraint(GRBModel& model, const std::vector<Node>& subproblem_nodes, const std::vector<uv>& pair_keys,
                      std::map<uv, GRBVar>& x_map ,std::map<ouv, GRBVar>& y_map)
{
    for (auto o : subproblem_nodes)
    {
        auto o_id = Graph::id(o);
        // For each x_{i,j}: using pair_keys we mantain the v > u restriction
        for (auto pair : pair_keys)
        {
            int u_id = std::get<0>(pair);
            int v_id = std::get<1>(pair);
            ouv l(o_id,u_id,v_id);
            ouv r(o_id,v_id,u_id);
            model.addConstr(y_map[l] + y_map[r], GRB_LESS_EQUAL, x_map[pair]);
        }
    }
}

// Avoid cicle constraint
// Similar to fourth constraint
void sixth_constraint(GRBModel& model, const std::vector<Node>& subproblem_nodes,
                      const std::vector<uv>& pair_keys, std::map<uv, GRBVar>& x_map)
{
    GRBLinExpr cicle_constr_linexp(0);
    // constraint: sum(x) <= n-1
    for (auto pair : pair_keys)
    {
        cicle_constr_linexp += x_map[pair];
    }
    model.addConstr(cicle_constr_linexp, GRB_EQUAL, (subproblem_nodes.size()-1));
}


Node findCentroidOfTree(Node n, SubTree& subtree, std::vector<Node> subproblem_nodes)
{
    Node min_node = n;
    int min_dist = INT32_MAX;
    for (auto u : subproblem_nodes)
    {
        int distance = 0;
        lemon::Dfs<SubTree> dfs(subtree);
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
    for (SubTree::IncEdgeIt it(subtree,u); it != lemon::INVALID; ++it)
    {
        Node v = subtree.oppositeNode(u, it);
        if (v == from) continue;
        countNodesFromEdge(v, u, subtree, nodes_from_edge);
        if (curr_edge == lemon::INVALID) continue;
        nodes_from_edge[curr_edge] += nodes_from_edge[it];
    }
}

void clusterFromTreeByEdge(Node n, Edge last_edge, SubTree& subtree, std::vector<Node>& cluster)
{
    cluster.push_back(n);
    for (SubTree::IncEdgeIt it(subtree,n); it != lemon::INVALID; ++it)
    {
        if (it == last_edge) continue;
        Node u = subtree.oppositeNode(n,it);
        clusterFromTreeByEdge(u, it, subtree, cluster);
    }
}


// Split tree in two if contains a perfect centroid (tree_size/2)
bool splitInEdge(int tree_size, SubTree& subtree, SubTree::EdgeMap<int>& nodes_from_edge, int i, int j,
                 Solution* solution)
{
    // Iterate over edges to find a edge centroid
    for (SubTree::EdgeIt it(subtree); it != lemon::INVALID; ++it)
    {
        if (nodes_from_edge[it] == ceil(tree_size / 2))
        {
            // Divide in this edge
            Node u = subtree.u(it);
            Node v = subtree.v(it);
            std::vector<Node> cluster_u{};
            std::vector<Node> cluster_v{};
            clusterFromTreeByEdge(u, it, subtree, cluster_u);
            clusterFromTreeByEdge(v, it, subtree, cluster_v);
            solution->clusters[i] = cluster_u;
            solution->clusters[j] = cluster_v;
            return true;
        }
    }
    return false;
}

void knapsackApproach(Node centroid, SubTree& subtree, SubTree::EdgeMap<int>& nodes_from_edge,
                      int tree_size, int i, int j, Instance& instance, Solution* solution, GRBEnv& env)
{
    GRBModel model_kp {GRBModel(env)};
    std::vector<GRBVar> vars;
    GRBLinExpr objective_expr(0);
    GRBLinExpr constraint_expr(0);

    for (SubTree::IncEdgeIt it(subtree,centroid); it != lemon::INVALID; ++it)
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
        std::vector<Node> cluster_one{centroid};
        std::vector<Node> cluster_two{centroid};
        for (SubTree::IncEdgeIt it(subtree,centroid); it != lemon::INVALID; ++it)
        {
            Node u = subtree.oppositeNode(centroid, it);
            if (u == lemon::INVALID) continue;

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
        solution->clusters[i] = cluster_one;
        solution->clusters[j] = cluster_two;
        // cout << "new1 = " << cluster_one.size() << "; new2 = " << cluster_two.size() << endl;
        // if (abs(((int)cluster_one.size()) - ((int)cluster_two.size())) > 2 )
        // {
        //     cout << "centroid: " << graph.id(centroid) << "; cluster lengths: " << ceil(tree_size/2) << endl;
        //     printEdgesSubTree(subtree);
        // }
        // cout << "==========" << endl;
    }
}

void splitTree(SubTree& subtree, std::vector<Node>& subproblem_nodes, int i, int j,
               Instance& instance, Solution* solution, GRBEnv& env, int debug=false)
{
    // first of all, find a centroid
    Node centroid = findCentroidOfTree(subproblem_nodes[0], subtree, subproblem_nodes);

    // After, count nodes hanged in each edge
    SubTree::EdgeMap<int> nodes_from_edge(subtree, 1);
    countNodesFromEdge(centroid, lemon::INVALID, subtree, nodes_from_edge);

    // After that, try to split the tree in an edge, because it's simple
    bool divided = splitInEdge(subproblem_nodes.size(), subtree, nodes_from_edge, i, j, solution);
    if (divided)
    {
        if (debug)
            std::cout << "divided in edge" << std::endl;
        return;
    }

    // Also, try to split in Node
    // divided = splitInNode(centroid, subtree, i, j);
    if (divided)
    {
        if (debug)
            std::cout << "divided in centroid" << std::endl;
        return;
    }

    // If nothing works, use a knapsack approach
    knapsackApproach(centroid, subtree, nodes_from_edge, subproblem_nodes.size(), i, j, instance, solution, env);
    if (debug)
        std::cout << "divided using knapsack approach" << std::endl;

}

double solveSubproblem(std::vector<Node> subproblem_nodes, bool* was_modified, int i, int j,
                       Instance& instance, Solution* solution, GRBEnv env, bool debug = false)
{
    // Cria nova tabela de requisitos copiando os valores originais
    auto subproblem_requirements {generateSubproblemsReq(subproblem_nodes, instance, solution)};

    double init_value = calculateSubproblemObjective(subproblem_nodes, subproblem_requirements, instance, solution);
    try
    {
        GRBModel model {GRBModel(env)};

        std::vector<uv> pair_keys {};
        std::map<uv, GRBVar> x_map;

        // f_{o,u,v} means: the flow of o going from u to v
        std::vector<ouv> triple_keys {};
        std::map<ouv, GRBVar> f_map;
        std::map<ouv, GRBVar> y_map;

        initVars(model, subproblem_nodes, pair_keys, triple_keys, x_map, f_map, y_map, instance);

        // First of all, the objective
        objective(model, subproblem_nodes, f_map, instance);

        // The origin sends the sum of all its requirements as initial flow
        first_constraint(model, subproblem_nodes, subproblem_requirements, f_map, instance);

        // Flow conservation
        second_constraint(model, subproblem_nodes, subproblem_requirements, f_map, instance);

        // f(o,u,v) must b zero if arc is not used
        third_constraint(model, subproblem_nodes, subproblem_requirements, f_map, y_map, instance);

        // Tree constraint for y
        fourth_constraint(model, subproblem_nodes, y_map, instance);

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

            if (!verifyTree(instance, solution, subproblem_nodes[0]))
            {
                if (debug)
                    perror("new solution not valid!\n");
                return 0;
            }
            if (((int)solver_result) < ((int) init_value)) // Update tree only if we got a better solution
            {
                if (debug)
                {
                    std::cout << "optimized!" << std::endl;
                    std::cout << "Subproblem result: " << solver_result;
                    std::cout << "; Original result: " << init_value << std::endl;
                }
                for (auto pair : pair_keys)
                {
                    auto x = x_map[pair];
                    auto x_name  = x.get(GRB_StringAttr_VarName);
                    auto x_value = (bool) x.get(GRB_DoubleAttr_X);
                    auto u = instance.graph.nodeFromId(std::get<0>(pair));
                    auto v = instance.graph.nodeFromId(std::get<1>(pair));
                    auto e = findEdge(instance.graph, u, v);
                    if (e != lemon::INVALID)
                    {
                        if (solution->edge_in_tree->operator[](e) != (bool) x_value)
                        {
                            *was_modified = true;
                        }
                        solution->edge_in_tree->operator[](e) = x_value;
                    }
                }

                solution->tree = new Tree(instance.graph, *solution->edge_in_tree);
                TreeNodeMapBool in_subtree(*solution->tree, false);

                for (auto n : subproblem_nodes)
                {
                    in_subtree[n] = true;
                }
                SubTree subtree(*solution->tree, in_subtree);
                auto calcted_value = 1;// calculateSubproblemObjective(subproblem_nodes, subproblem_requirements);
                // printEdgesSubTree(subtree);
                splitTree(subtree, subproblem_nodes, i, j, instance, solution, env);

                // if (solver_result != calcted_value)
                // {
                std::cout << "initial: " << init_value << std::endl;
                std::cout << "c1 = " << i << "; c2 = " << j << std::endl;
                std::cout << "solver: " << solver_result << "\tcalculated: " << calcted_value << std::endl;
                printClusters(solution);
                // printEdgesTree();
                std::cout << "----------------------------------------------------------" << std::endl;
                // }

                return solver_result - init_value;
            }
            else
            {
                if (debug)
                {
                    std::cout << "no optimized" << std::endl
                         << "Subproblem result: " << solver_result
                         << "; Original result: " << init_value << std::endl;
                }
                return 0;
            }
        }
        else {
            if (debug)
                std::cout << "fail" << std::endl;
            return 0;
        }

    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    return 0;
}

#endif //MATH_HEURISTIC_HEURISTIC_H
