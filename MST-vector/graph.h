#pragma once
#include <fstream>
#include <vector>
#include <unordered_map>
#include "hash.hpp"
#include "utility.hpp"
#include "time.h"
#include <iostream>

class Graph{
    public:
        int nLabels, n, m, max_core, maxDeg;
        std::vector<CuckooHash> nbrMap;
        std::vector<int> degree;
        std::vector<std::vector<int>> nbr;
        std::unordered_map<int,int> labelId;
        std::vector<int> vLabel;
        std::vector<int> pstart,edges,pend;
        std::vector<int> order, ordered, core;

        Graph(const std::string& dataset) {
            clear();
            loadFromFile(dataset);
        }

        Graph() {
            clear();
        }

        Graph(const Graph& graph) = default;
        Graph(Graph&& graph) = default;
        Graph& operator=(const Graph& graph) = default;
        Graph& operator=(Graph&& graph) = default;

        ~Graph() {
            clear();
        }

        void loadFromFile(const std::string& filename){
            printf("# Start reading graph %s\n", filename.c_str());
            FILE *f = Utility::open_file(filename.c_str(), "r");
            fscanf(f, "%d", &nLabels);
            int numOfEdges;
            fscanf(f, "%d", &numOfEdges);

            for(int i = 0;i<numOfEdges;++i){
                int u, v;
                fscanf(f, "%d", &u);
                fscanf(f, "%d", &v);
                addEdgeWithLabel(u, v);
            }
            
            fclose(f);
            printf("There are %d nodes and %d edges \n", n, m);

			toAdjArr();  
        }


		void toAdjArr(){
            pstart.reserve(n+1);
            pend.reserve(n);
            edges.reserve(2*m);
            int numEdges = 0;
            for(int u = 0;u<n;++u){
                pstart[u] = u==0 ? 0 : pend[u-1] + 1;
                for(auto v: nbr[u]){
                    edges[numEdges++] = v;
                }
                pend[u] = pstart[u]+degree[u]-1;
            }
            pstart[n] = pend[n-1] + 1;
		}

        void clear() {
            m = 0;
            maxDeg = nLabels = n = 0;
            degree.clear();
            nbr.clear();
            nbrMap.clear();
            vLabel.clear();
            labelId.clear();
            pstart.clear();
            edges.clear();
            pend.clear();
            order.clear();
            ordered.clear();
            core.clear();
        }

        void removeEdge(int u, int vIdx){
            std::swap(edges[vIdx], edges[pend[u]]);
            --pend[u];
        }

        int getIndex(int label) {
            if (labelId.count(label)) return labelId[label];
            vLabel.push_back(label);
            degree.push_back(0);
            nbr.push_back(std::vector<int>());
            nbrMap.push_back(CuckooHash());
            return labelId[label] = n++;
        }

        void addEdge(int u, int v) {
            if (nbrMap[u].find(v))return;
            ++m;
            nbr[u].push_back(v);
            nbr[v].push_back(u);
            nbrMap[u].insert(v);
            nbrMap[v].insert(u);
            ++degree[u];
            ++degree[v];
            maxDeg = std::max(maxDeg, degree[u]);
            maxDeg = std::max(maxDeg, degree[v]);
        }

        void addEdgeWithLabel(int ulabel, int vlabel) {
            nLabels = std::max(nLabels, ulabel + 1);
            nLabels = std::max(nLabels, vlabel + 1);
            addEdge(getIndex(ulabel), getIndex(vlabel));
        }

        void reorganizeAdjacencyList(){
            for(int i = 0;i<n;i++){
                pend[i] = pstart[i] - 1;
                for(int j = pstart[i];j<pstart[i+1];j++){
                    if(order[edges[j]]>order[i]){
                        edges[++pend[i]] = edges[j];
                    }
                }
            }
        }

        void print(){
            std::cerr << "print graph: " << std::endl;
            for(int i = 0;i<n;i++){
                for(int j = pstart[i];j<=pend[i];j++){
                    std::cerr << i << " " << edges[j] << std::endl;
                }
            }
        }
};