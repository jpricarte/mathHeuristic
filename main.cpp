#include <iostream>
#include <fstream>
#include <chrono>
#include "instance.h"
#include "heuristic.h"

using namespace std;

GRBEnv env(false);

int main (int argc, char* argv[]) {

    if (argc != 4)
    {
        perror("usage: ./heuristic inputFile clusterSize maxIter \n");
        return 1;
    }
    env.set("OutputFlag", "0");
    env.start();

    string input_file = argv[1];
    int cluster_size = stoi(argv[2]);
    int max_iteration = stoi(argv[3]);
    Instance instance{};
    readInstance(input_file, instance);
    auto start = chrono::high_resolution_clock::now();
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
    printEdgesTree(instance,solution);

    int iterNum = 0;
    vector<int> modified_clusters = {};

    for (int i = 0; i < solution->clusters.size(); i++)
    {
        cout << modified_clusters.size() << endl;
        for (int j = i+1; j < solution->clusters.size(); j++) {
            auto some_cluster = selectTwoClusters(solution, i, j);
            if (some_cluster.empty())
                continue;
            iterNum++;
            bool was_modified = false;
            auto diff = solveSubproblem(some_cluster, &was_modified, i, j, instance, solution, env);
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
        }
    }

    while (!modified_clusters.empty()) {
        int i = modified_clusters.front();
        modified_clusters.erase(modified_clusters.begin());
        for (int j = 0; j < solution->clusters.size(); j++) {
            if (i == j) continue;

            auto some_cluster = selectTwoClusters(solution, i, j);
            if (some_cluster.empty())
                continue;
            iterNum++;

            bool was_modified = false;

            // return the diference between the original value and the new solution
            auto diff = solveSubproblem(some_cluster, &was_modified, i, j, instance, solution, env);


            auto new_objective = objective + diff;
            if (new_objective < objective)


                if (was_modified) {
                    if (find(modified_clusters.begin(), modified_clusters.end(), i) == modified_clusters.end()) {
                        modified_clusters.push_back(i);
                    }
                    if (find(modified_clusters.begin(), modified_clusters.end(), j) == modified_clusters.end()) {
                        modified_clusters.push_back(j);
                    }
                }
        }
    }

    objective = calculateObjective(instance, solution);
    std::cout << "final objective: " << objective << std::endl;
    auto stop = chrono::high_resolution_clock::now();
    auto exec_time = chrono::duration_cast<chrono::seconds>(stop - start);
    auto name_index = input_file.find_last_of("/");
    auto ext_index = input_file.find_last_of(".");
    auto stats_file = "./output/" +
                      input_file.substr(name_index, ext_index) +
                      "_stats.csv";
    std::ofstream outfile;
    outfile.open(stats_file, std::ios_base::app); // append instead of overwrite
    outfile << cluster_size << "," << solution->clusters.size() << iterNum << ","
            << exec_time.count() << "," << objective << endl;
    outfile.close();
    return 0;
}
