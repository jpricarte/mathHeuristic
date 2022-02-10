#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>

using namespace std;
template <typename T>
class Graph{
public:
    Graph(int nEdge);
    T get(int u, int v);
    void set(int u, int v, T val);
    vector<T> graphList;
    int size;
};

template <typename T>
Graph<T>::Graph(int nEdge) : size(nEdge) graphList(vector<T>(nEdge)) {}

template <typename T>
T Graph<T>::get(int u, int v) {
    if (v > u) swap(u,v);
    return this->graphList[floor(((u+1)*u)/2)+v];
}

template <typename T>
void Graph<T>::set(int u, int v, T val) {
    if (v > u) swap(u,v);
    this->graphList[floor(((u+1)*u)/2)+v] = val;
}


typedef struct NodePair {
    double req;
    double length;
    bool inTree;

    NodePair() : req(0.0), length(-1), inTree(false) {}
    NodePair(double _req, double _length, bool _inTree) : req(_req), length(_length), inTree(_inTree) {}
} Edge;


Graph<Edge*> graph(1);
int n, m;

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
    graph = Graph<Edge*>(int((n+1)*n/2));
    for (uint32_t i=0; i<graph.graphList.size(); i++) {
        graph.graphList[i] = new Edge();
    }

    // Adding edges and length
    for (auto i=0; i<m; i++)
    {
        instanceFile >> e1 >> e2 >> e3;
        u = stoi(e1);
        v = stoi(e2);
        value = stod(e3);
        graph.get(u,v)->length = value;
    }
    // Adding requirements
    for (auto i=0; i<n; i++)
    {
        for (auto j=i; j<n; j++)
        {
            instanceFile >> e1;
            value = stod(e1);
            graph.get(j,i)->req = value;
        }
    }
    instanceFile.close();
}

int main(int argc, char* argv[]) {
    if (argc != 2)
    {
        cout << "usage: ./heuristic inputFile" << endl;
    }

    // Inicialization
    // Read instance
    readInstance(argv[1]);

}