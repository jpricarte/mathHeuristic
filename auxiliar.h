//
// Created by jpricarte on 06/09/22.
//

#ifndef MATH_HEURISTIC_AUXILIAR_H
#define MATH_HEURISTIC_AUXILIAR_H

#include "instance.h"

void printEdgesTree(const Instance& instance, const Solution* solution)
{
    for (auto edge : instance.edges)
    {
        if (solution->edge_in_tree->operator[](edge))
        {
            std::cout << Graph::id(instance.graph.u(edge))
                      << " " << Graph::id(instance.graph.v(edge)) << std::endl;
        }
    }
}

bool containsCycle(const Node n, const Node up, TreeNodeMapBool& visited, const Solution* solution)
{
    bool has_cycle = false;
    visited[n] = true;
    for (Tree::IncEdgeIt e(*(solution->tree), n); e != lemon::INVALID; ++e)
    {
        Node v = solution->tree->oppositeNode(n,e);
        if (v != up)
        {
            if (visited[v])
            {
                return true;
            }
            has_cycle = containsCycle(v, n, visited, solution);
            // If some subgraph return a cycle, break the algorithm and return true
            if (has_cycle) return has_cycle;
        }
    }
    return false;
}

void printClusters(Solution* solution)
{
    for (const auto& cluster : solution->clusters)
    {
        std::cout << "{ ";
        for (auto node : cluster)
        {
            std::cout << solution->tree->id(node) << " ";
        }
        std::cout << "}" << std::endl;
    }
}

bool verifyTree(const Instance& instance, const Solution* solution, const Node root)
{
    TreeNodeMapBool visited(*(solution->tree), false);
    bool contains_cycle = containsCycle(root, lemon::INVALID, visited, solution);
    if (contains_cycle)
    {
        std::cout << "graph contains cycle" << std::endl;
        return false;
    }

    bool connected = true;
    for (auto node : instance.nodes)
    {
        if (!visited[node])
        {
            connected = false;
            std::cerr << "graph contains more than one connected component" << std::endl;
            break;
        }
    }
    // Is a tree if not contains cycle and is connected
    return connected;
}

std::vector<Node> selectTwoClusters(Solution* solution, int first, int second)
{
    std::vector<Node> final_cluster = {};

    for (auto node_one : solution->clusters[first])
    {
        for (auto node_two : solution->clusters[second])
        {
            // If clusters are connected, create a merged cluster
            if (node_one == node_two || findEdge(*solution->tree, node_one, node_two) != lemon::INVALID)
            {
                for (auto node : solution->clusters[first])
                {
                    final_cluster.push_back(node);
                }

                for (auto node : solution->clusters[second])
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

// Aux function for subproblem requirements
/*
    base-index: indice no vetor do subproblema da aresta
a árvore que será adicionada
    current: vértice que está sendo analizado nesse momento
    previous: último vértice analisado, apenas para evitar volta
    outros: vetores comuns
*/
void addRequirements(int base_index, Node current , Node previous,
                     std::vector<Node> subproblem_nodes,
                     std::vector<std::vector<double>> *subproblem_requirements,
                     Instance& instance, Solution* solution)
{
    auto current_id = Graph::id(current);

    // Para cada vértice da árvore
    for (auto i=0; i < (int) subproblem_nodes.size(); i++)
    {
        // Inicialmente, o requerimento entre o vértice da árvore que será sobrecarregado (base_index)
        // e os outros vértices receberá uma soma do requerimento do novo vértice (current)
        auto node_id = Graph::id(subproblem_nodes[i]);
        if (i > base_index)
        {
            (*subproblem_requirements)[base_index][i] += instance.requirements[node_id][current_id];
        }
        else
        {
            (*subproblem_requirements)[i][base_index] += instance.requirements[node_id][current_id];
        }
    }

    // Recursivamente, vai para os outros vértices ligados a esse que ainda não foram acessados
    // Por ser uma árvore, não precisamos de controle de acessados
    for (auto e = Tree::IncEdgeIt(*solution->tree,current); e != lemon::INVALID; ++e)
    {
        auto node = solution->tree->oppositeNode(current, e);
        if (node == previous)
        {
            continue;
        }
        addRequirements(base_index, node, current, subproblem_nodes,
                        subproblem_requirements, instance, solution);
    }
}

std::vector<std::vector<double>> generateSubproblemsReq(const std::vector<Node>& subproblem_nodes,
                                              Instance& instance, Solution* solution)
{
    std::vector<std::vector<double>> subproblem_requirements = {};
    for (auto i=0; i < (int) subproblem_nodes.size(); i++)
    {
        subproblem_requirements.emplace_back();

        // Cria matriz auxiliar de requerimentos (triangular superior)
        for (auto j=0; j < (int) subproblem_nodes.size(); j++)
        {
            if (j > i)
            {
                auto u_id = Graph::id(subproblem_nodes[i]);
                auto v_id = Graph::id(subproblem_nodes[j]);
                subproblem_requirements[i].push_back(instance.requirements[u_id][v_id]);
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
        for (auto e = Tree::IncEdgeIt(*solution->tree,u); e != lemon::INVALID; ++e)
        {
            auto v = solution->tree->oppositeNode(u, e);
            // se o vertice oposto não estiver no subproblema,
            // vai somando todos os requisitos dos vértices agregados para ele
            if (find(subproblem_nodes.begin(), subproblem_nodes.end(), v) == subproblem_nodes.end())
            {
                addRequirements(i, v, u, subproblem_nodes,
                                &subproblem_requirements, instance, solution);
            }
        }
    }

    for (auto i=0; i < (int) subproblem_nodes.size(); i++)
    {
        subproblem_requirements[i][i] = 0;
    }

    return subproblem_requirements;
}

#endif //MATH_HEURISTIC_AUXILIAR_H
