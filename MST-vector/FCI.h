#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "graph.h"
#include <cassert>
#include <algorithm>
#include "heap.h"
#include "disjointSet.h"
#include "tree.h"
#include "tree1.h"
#include <list>
#include <queue>
#include "Timer.h"
#include <cmath>

class Cluster{
public:
    Cluster(int rep_): rep(rep_){}
    int rep;
    std::vector<int> cores, nonCores;
};

class EdgeNode{
public:
    Edge edge;
    double bestDistance;
};

class VertexInfo{
public:
    int id;
    double distance;
    bool operator<(const VertexInfo& other) const {
        return distance > other.distance;
    }
};

class FCI{
public:
    int mu;
    Graph g;
    int n1, m;
    std::vector<VertexInfo> vertices;
    std::vector<double> distance_threshold;
    std::vector<int> count, two_hop_neighbors;
    std::vector<std::vector<Edge>> FCF;
    std::vector<int> rFCF;
    std::vector<Edge> NCI;
    std::vector<Tree> FCFTree;
    std::vector<Tree1> FCFTree1;
    std::vector<Edge> NCIEdges;
    std::vector<Edge> FCFEdges;
    std::vector<std::pair<int, std::vector<double>>> data;
    

    void initialise(){
        vertices.clear();
        distance_threshold.clear();
        distance_threshold.assign(n1,0);
        rFCF.assign(n1, -1);
        count.clear();
        count.assign(n1, 0);
        two_hop_neighbors.clear();
        FCF.clear();
        NCI.clear();
        FCFTree.clear();
        FCFTree1.clear();
        FCFEdges.clear();
    }

        
    std::vector<Edge> read_edges(const std::string &file_name, const double &eps) {
        std::vector<Edge> edges;
        std::ifstream file(file_name);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file " << file_name << std::endl;
            return edges;
        }
    
        int u, v;
        double weight;
        while (file >> u >> v >> weight) {
            if (weight > eps) break;
            edges.push_back({u, v, weight});
        }
    
        file.close();
        return edges;
    }

    void read_vectors(const std::string &file_name) {
        std::ifstream fin(file_name);
        if (!fin.is_open()) {
            std::cerr << "Failed to open file: " << file_name << std::endl;
            return;
        }
    
        std::string line;
    
        if (std::getline(fin, line)) {
            std::stringstream ss(line);
            ss >> n1;
        }
    
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
    
            std::stringstream ss(line);
            int index;
            char colon;
    
            ss >> index >> colon;
            if (colon != ':') continue;
    
            std::vector<double> vec;
            double val;
            while (ss >> val) {
                vec.push_back(val);
            }
    
            if (!vec.empty()) data.push_back({index, std::move(vec)});
        }
    
        distance_threshold.clear();
        distance_threshold.assign(n1, 0);
        rFCF.assign(n1, -1);
        count.clear();
        count.assign(n1, 0);
    }

    void bulid_core_index(){
        distance_threshold.resize(data.size());
        vertices.clear();
    
        std::priority_queue<double, std::vector<double>, std::less<double>> pq; 

        for (int u = 0; u < data.size(); ++u) {
            const auto &vec_u = data[u].second;
            for (int w = 0; w < data.size(); ++w) {
                if (w == u) continue;
                const auto &vec_w = data[w].second;
    
                double dist = 0.0;
                for (int k = 0; k < vec_u.size(); ++k) {
                    double diff = vec_u[k] - vec_w[k];
                    dist += diff * diff;
                }
                dist = std::sqrt(dist);
                if (pq.size() < (int)(mu - 1)) {
                    pq.push(dist);
                } else if (dist < pq.top()) {
                    pq.pop();
                    pq.push(dist);
                }
            }
    
            double dist_threshold = pq.size() >= (int)(mu - 1) ? pq.top() : 0;
            distance_threshold[u] = dist_threshold;
            vertices.push_back({data[u].first, dist_threshold});
            while (!pq.empty()) pq.pop();
        }
        std::sort(vertices.begin(), vertices.end());
    }

    void build_spanning_tree(){
        n1 = data.size();
        Heap heap = Heap(n1);
        std::vector<bool> visited(n1, false);

        std::vector<Edge> currEdges;
        for(auto vertex: vertices){
            int u = vertex.id;
            if(visited[u]) continue;

            double lowest_distance = 9999999;
            int closest_core = -1;

            for(auto vec:data){
                double dist = 0.0;
                int w = vec.first;
                for (int k = 0; k < vec.second.size(); ++k) {
                    double diff = data[u].second[k] - vec.second[k];
                    dist += diff * diff;
                }
                dist = std::sqrt(dist);
                double weight = std::max(dist, distance_threshold[w]);
                if((distance_threshold[w] < distance_threshold[u]) && (lowest_distance >= weight)){
                    lowest_distance = weight;
                    closest_core = w;
                }
                if(!visited[w]){
                    heap.addToHeap({w, u, std::max(weight, distance_threshold[u])});
                }
            }
            if(closest_core!=-1 && lowest_distance < distance_threshold[u]) NCI.push_back({closest_core, u,lowest_distance});
            visited[u] = true;
            rFCF[u] = FCF.size();
            

            while(heap.size()){
                auto cur = heap.pop();
                int u = cur.id;
                if(visited[u]) continue;
                visited[u] = true;
                rFCF[u] = FCF.size();
                currEdges.push_back({u,cur.id2,cur.weight});

                double lowest_distance = 999999;
                int closest_core = -1;

                for(auto vec:data){
                    double dist = 0.0;
                    int w = vec.first;
                    for (int k = 0; k < vec.second.size(); ++k) {
                        double diff = data[u].second[k] - vec.second[k];
                        dist += diff * diff;
                    }
                    dist = std::sqrt(dist);

                    auto weight = std::max(dist, distance_threshold[w]);
                    if((distance_threshold[w] < distance_threshold[u]) && (lowest_distance >= weight)){
                        lowest_distance = weight;
                        closest_core = w;
                    }
                    if(visited[w])continue;
                    if(!heap.inside(w)){
                        heap.addToHeap({w, u, std::max(weight, distance_threshold[u])});
                    }else{
                        heap.increasePriority({w, u, std::max(weight, distance_threshold[u])});
                    }
                }
               
                if(closest_core!=-1  && lowest_distance<distance_threshold[u]) NCI.push_back({closest_core, u,lowest_distance}); //NCI {core, noncore}
            }
            FCF.emplace_back(std::move(currEdges));
            currEdges.clear();
        }
        std::sort(NCI.begin(), NCI.end());
    }

    void build_NCI_density(){
        for(auto edges: FCF) FCFTree.push_back(Tree(edges));
        long long allEdges = 0, restEdges = 0;
        std::vector<Edge> edges, compulsoryEdges;
        
        for(auto vec:data){
            int u = vec.first;
            for(auto vec2:data){
                double dist = 0.0;
                int w = vec2.first;
                for (int k = 0; k < vec2.second.size(); ++k) {
                    double diff = vec.second[k] - vec2.second[k];
                    dist += diff * diff;
                }
                dist = std::sqrt(dist);

                auto weight = std::max(dist, distance_threshold[w]);
                if(weight < distance_threshold[u]){
                    edges.push_back({w,u,weight});             // NCIEdges {core, noncore}
                }
            }
            
            sort(edges.begin(), edges.end());
            int treeIdx = rFCF[u];
            auto &tree = FCFTree[treeIdx];
            for(int i = 0;i<edges.size();i++){
                int v = edges[i].id;
                bool pruned = false;
                for(int j = 0;j<compulsoryEdges.size();j++){
                    int w = compulsoryEdges[j].id;
                    auto maxWeight = tree.queryPathMin(tree.getIndex(v), tree.getIndex(w));
                    if(maxWeight <= edges[i].weight){
                        pruned = true;
                        break;
                    }
                }
                if(!pruned) {
                    compulsoryEdges.push_back(edges[i]);
                }
            }
            allEdges += edges.size(), restEdges += compulsoryEdges.size();
            NCIEdges.insert(NCIEdges.end(),
                std::make_move_iterator(compulsoryEdges.begin()),
                std::make_move_iterator(compulsoryEdges.end()));
            compulsoryEdges.clear();
            edges.clear();
        }
        std::sort(NCIEdges.begin(), NCIEdges.end());
        std::cout << restEdges << std::endl;
        std::cout << restEdges << " " << allEdges << std::endl;
        std::cout << NCI.size() << std::endl;
    }

    std::vector<Cluster> exactClustering(const std::string &directory_name, int mu, const double& eps){
        FCFEdges = read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"newFCI_sorted.txt", eps);
        NCIEdges = read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"newNCI_sorted.txt", eps);
        std::ofstream log_file("clustering_log.txt", std::ios::app);
        Timer t;
        std::vector<Cluster> result;
        DisJointSet clusters = DisJointSet(n1);
        std::vector<bool> visited = std::vector<bool>(n1, false);
        std::vector<int> vertex2Cluster = std::vector<int>(n1, -1);
        for(auto &edge: FCFEdges){
            int u = edge.id, v = edge.id2;
            clusters.my_union(u, v);
        }

        for(auto &edge: FCFEdges){
            int v = edge.id, u = edge.id2;
            int pv = clusters.find_root(v);
            if(vertex2Cluster[pv]==-1){
                vertex2Cluster[pv] = result.size();
                result.push_back(Cluster(pv));
            }
            if(!visited[v]){
                result[vertex2Cluster[pv]].cores.push_back(v);
                visited[v] = true;
            }
            if(!visited[u]){
                result[vertex2Cluster[pv]].cores.push_back(u);
                visited[u] = true;
            }
        }

        for(auto &edge: NCIEdges){
            int v = edge.id, u = edge.id2; // core noncore
            int pv = clusters.find_root(v),  pu = clusters.find_root(u);
            if(vertex2Cluster[pv]==-1){
                vertex2Cluster[pv] = result.size();
                result.push_back(Cluster(pv));   
                result[vertex2Cluster[pv]].cores.push_back(v);            
            }
             if(vertex2Cluster[pu]==-1){
                int idx = vertex2Cluster[pv];
                result[idx].nonCores.push_back(u);
            }  
        }

        auto elapsed = t.elapsed();
            std::string msg = "Exact | Dataset: " + directory_name + ", mu: " + std::to_string(mu) +
                            ", eps: " + std::to_string(eps) +
                            ", time(ms): " + Utility::integer_to_string(elapsed);
            std::cout << msg << std::endl;
            log_file << msg << std::endl;
            log_file.close();
        return result;
    }

    std::vector<Cluster> densityClustering(const std::string &directory_name, int mu, const double& eps){
        FCFEdges = std::move(read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"newFCI_sorted.txt", eps));
        NCIEdges = std::move(read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"newDNCI_sorted.txt", eps));
        std::ofstream log_file("clustering_log.txt", std::ios::app);
        Timer t;
        std::vector<Cluster> result;
        DisJointSet clusters = DisJointSet(n1);
        std::vector<int> vertex2Cluster = std::vector<int>(n1, -1);
        std::vector<bool> visited = std::vector<bool>(n1, false);
        for(auto &edge: FCFEdges){
            int u = edge.id, v = edge.id2;
            clusters.my_union(u, v);
        }
        for(auto &edge: FCFEdges){
            int v = edge.id, u = edge.id2;
            int pv = clusters.find_root(v);
            if(vertex2Cluster[pv]==-1){
                vertex2Cluster[pv] = result.size();
                result.push_back(Cluster(pv));
            }
            if(!visited[v]){
                result[vertex2Cluster[pv]].cores.push_back(v);
                visited[v] = true;
            }
            if(!visited[u]){
                result[vertex2Cluster[pv]].cores.push_back(u);
                visited[u] = true;
            }
        }

        for(auto &edge: NCIEdges){
            int v = edge.id, u = edge.id2; // core noncore
            int pv = clusters.find_root(v),  pu = clusters.find_root(u);
            if(vertex2Cluster[pv]==-1){
                vertex2Cluster[pv] = result.size();
                result.push_back(Cluster(pv));   
                result[vertex2Cluster[pv]].cores.push_back(v);            
            }
             if(vertex2Cluster[pu]==-1){
                int idx = vertex2Cluster[pv];
                result[idx].nonCores.push_back(u);
            }  
        }
        for(auto &cluster: result){
            std::vector<int> tmp;
            for(auto v: cluster.nonCores){
                if(!visited[v]){
                    tmp.push_back(v);
                    visited[v] = true;
                }
            }
            for(auto v: tmp)visited[v] = false;
            cluster.nonCores = std::move(tmp);
        }

        auto elapsed = t.elapsed();
            std::string msg = "Density | Dataset: " + directory_name + ", mu: " + std::to_string(mu) +
                            ", eps: " + std::to_string(eps) +
                            ", time(ms): " + Utility::integer_to_string(elapsed);
            std::cout << msg << std::endl;
            log_file << msg << std::endl;
            log_file.close();
        return result;
    }
};