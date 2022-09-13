#include <iostream>
#include <fstream>
#include "instance.h"
#include "heuristic.h"

using namespace std;

GRBEnv env(false);

void readInstance(const string& filename, Instance& instance);

int main (int argc, char* argv[]) {
    if (argc != 4)
    {
        perror("usage: ./heuristic inputFile clusterSize maxIter \n");
        return 1;
    }
    string input_file = argv[1];
    int cluster_size = stoi(argv[2]);
    int max_iteration = stoi(argv[3]);
    Instance instance{};
    readInstance(input_file, instance);
    auto root = chooseRoot(instance, cluster_size);
    auto solution = generateInitialSolutionDijkstra(instance, root, cluster_size);
    if (!verifyTree(instance, solution, root))
    {
        perror("Initial solution not valid!\n");
        printEdgesTree(instance, solution);
        return 1;
    }

    TreeNodeMapBool in_some_cluster{*solution->tree, false};
    bool could_divide = divideTree(root, lemon::INVALID, in_some_cluster, solution, cluster_size);
    if (!could_divide) {
        perror("Error in division");
        return 2;
    }
    if (solution->clusters.back().empty()) {
        solution->clusters.pop_back();
    }

    cout << "clusterized" << endl;
    cout << "clusters created: " << solution->clusters.size() << endl;
    printClusters(solution);
    double objective = calculateObjective(instance, solution);
    std::cout << "initial objective: " << objective << std::endl;

    int iterNum = 0;
    vector<int> modified_clusters = {};

    for (int i = 0; i < solution->clusters.size(); i++)
    {
        for (int j = i+1; j < solution->clusters.size(); j++) {
            auto some_cluster = selectTwoClusters(solution, i, j);
            if (some_cluster.empty())
                continue;
            iterNum++;
            bool was_modified = false;
            auto diff = solveSubproblem(some_cluster, &was_modified, i, j, instance, solution, env);
        }
    }

    // Currently it's just iterating over nodes, some approach will be used to converge in an optimum solution

    return 0;
}

void readInstance(const string& filename, Instance& instance)
{
    // Stream and temporary vars
    ifstream instanceFile(filename);
    if (!instanceFile.is_open())
    {
        cerr << "file not open" << endl;
        return;
    }
    string e1, e2, e3;
    unsigned int u, v;
    double value;
    // Get graph size (n for nodes and m for edges) and resizing everything
    instanceFile >> e1 >> e2;
    instance.n = stoi(e1);
    instance.m = stoi(e2);
    instance.graph.reserveNode(instance.n);
    instance.nodes.resize(instance.n);
    instance.requirements.resize(instance.n, vector<double>(instance.n));
    instance.graph.reserveEdge(instance.m);
    instance.edges.resize(instance.m);

    // Adding nodes
    for (auto i=0; i<instance.n; i++)
    {
        instance.nodes[i] = instance.graph.addNode();
    }

    // Adding edges and length
    for (auto i=0; i<instance.m; i++)
    {
        instanceFile >> e1 >> e2 >> e3;
        u = stoi(e1);
        v = stoi(e2);
        value = stod(e3);
        instance.edges[i] = instance.graph.addEdge(instance.nodes[u], instance.nodes[v]);
        instance.lengths[instance.edges[i]] = value;
    }
    // Adding requirements
    for (auto i=0; i<instance.n; i++)
    {
        instance.requirements[i][i] = 0;
        for (auto j=i+1; j<instance.n; j++)
        {
            instanceFile >> e1;
            value = stod(e1);
            instance.requirements[i][j] = value;
            instance.requirements[j][i] = value;
        }
    }
    instanceFile.close();

}
