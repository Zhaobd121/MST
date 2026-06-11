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

class Cluster{
public:
    Cluster(int rep_): rep(rep_){}
    int rep;
    std::vector<int> cores, nonCores;
};

class EdgeNode{
public:
    Edge edge;
    Fraction bestSimilarity;
};

class VertexInfo{
public:
    int id;
    Fraction similarity;
    bool operator<(const VertexInfo& other) const {
        return similarity > other.similarity;
    }
};

class FCI{
public:
    int mu;
    Graph g;
    int n1, n2, n, m;
    std::vector<VertexInfo> vertices;
    std::vector<Fraction> similarity_threshold;
    std::vector<int> count, two_hop_neighbors;
    std::vector<std::vector<Edge>> FCF;
    std::vector<int> rFCF;
    std::vector<Edge> NCI;
    std::vector<Tree> FCFTree;
    std::vector<Tree1> FCFTree1;
    std::vector<Edge> NCIEdges;
    std::vector<Edge> FCFEdges;
    

    void initialise(){
        vertices.clear();
        similarity_threshold.clear();
        similarity_threshold.assign(n1,{0,1});
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

    std::vector<Edge> read_edges(const std::string &file_name, const Fraction &eps){
        FILE* file = fopen(file_name.c_str(), "rb");
        
        int buffer[4];
        std::vector<Edge> edges;
        while (fread(buffer, sizeof(int), 4, file) == 4) {
            int u = buffer[0];
            int v = buffer[1];
            int num = buffer[2];
            int denom = buffer[3];
            if(Fraction{num, denom}<eps)break;
            edges.push_back({u,v,{num, denom}});
        }
        fclose(file);
        return edges;
    }


    void read_dataset(const std::string &file_name) {
        std::ifstream input_file(file_name, std::ios::in);
        input_file >> n1 >> n2 >> m;
        n = n1 + n2;
        std::cout << file_name << " : n1 = " << n1 << ", n2 = " << n2 << ", m = " << m << std::endl;
        std::string line;
        int id, value;
        for(int i = 0;i<n;i++) g.getIndex(i);
        while (getline(input_file, line)) {
            if (line.empty()) continue;

            std::stringstream ss(line);
            char colon;
            ss >> id >> colon;
            while (ss >> value) {
                assert(id-1==g.getIndex(id-1));
                assert(value + n1 -1==g.getIndex(value + n1 -1));
                g.addEdgeWithLabel(id-1, value+n1-1);
            }
        }
        input_file.close();
        assert(g.m==m);
        g.toAdjArr();
    }

    void bulid_core_index(int _mu){
        initialise();
        mu = _mu;
        std::vector<Fraction> two_hop_similarity;
        std::priority_queue<Fraction, std::vector<Fraction>, std::greater<Fraction> > pq;
        
        for(int u = 0; u<n1;u++){
            for(int i = g.pstart[u];i<g.pstart[u+1];i++){
                int v = g.edges[i];
                for(int j = g.pstart[v]; j<g.pstart[v+1];j++){
                    int w = g.edges[j];
                    if(w==u) continue;
                    if(count[w]==0) two_hop_neighbors.push_back(w);
                    count[w]++;
                }
            }

            for(auto w : two_hop_neighbors){
                Fraction sim = Fraction(count[w], g.degree[u] + g.degree[w] - count[w]);
                if (pq.size() < mu - 1) {
                    pq.push(sim);
                } else if (sim > pq.top()) {
                    pq.pop();
                    pq.push(sim);
                }
                count[w] = 0;
            }

            two_hop_neighbors.clear();
            // std::sort(two_hop_similarity.begin(),two_hop_similarity.end(),[](const Fraction &a, const Fraction &b){return a > b;});
            // auto sim = two_hop_similarity.size() >= mu - 1 ? two_hop_similarity[mu-2] : Fraction(0, 1);
            auto sim = pq.size() >= mu - 1 ? pq.top() : Fraction(0, 1);
            similarity_threshold[u] = sim;
            vertices.push_back({u, sim});
    
            while(pq.size())pq.pop();
            two_hop_similarity.clear();
        }
        std::sort(vertices.begin(),vertices.end());     
    }

    void build_spanning_tree(){
        Heap heap = Heap(n1);
        std::vector<bool> visited(n1, false);
        for(int i = 0;i<n1;i++)assert(count[i]==0);
        std::vector<Edge> currEdges;
        for(auto vertex: vertices){
            {
                int u = vertex.id;
                if(visited[u]) continue;
                for(int i = g.pstart[u]; i<g.pstart[u+1]; i++){
                    int v = g.edges[i];
                    for(int j = g.pstart[v]; j < g.pstart[v+1]; j++){
                        int w = g.edges[j];
                        if(w==u) continue;
                        if(count[w]==0) two_hop_neighbors.push_back(w);
                        count[w]++;
                    }
                }

                Fraction highest_similarity = {0, 1};
                int closest_core = -1;
                for(auto w: two_hop_neighbors){
                    auto sim = Fraction(count[w], g.degree[u] + g.degree[w] - count[w]);
                    auto weight = std::min(sim, similarity_threshold[w]);
                    if((similarity_threshold[w] > similarity_threshold[u]) && (highest_similarity <= weight)){
                        highest_similarity = weight;
                        closest_core = w;
                    }
                    if(!visited[w]){
                        heap.addToHeap({w, u, std::min(weight, similarity_threshold[u])});
                    }
                }
                for(auto w: two_hop_neighbors) count[w] = 0;
                two_hop_neighbors.clear();
                if(closest_core!=-1 && highest_similarity>similarity_threshold[u]) NCI.push_back({closest_core, u,highest_similarity});
                visited[u] = true;
                rFCF[u] = FCF.size();
            }
            

            while(heap.size()){
                auto cur = heap.pop();
                int u = cur.id;
                if(visited[u]) continue;
                visited[u] = true;
                rFCF[u] = FCF.size();
                currEdges.push_back({u,cur.id2,cur.weight});

                for(int i = g.pstart[u]; i < g.pstart[u+1]; i++){
                    int v = g.edges[i];
                    for(int j = g.pstart[v]; j < g.pstart[v+1]; j++){
                        int w = g.edges[j];
                        if(w == u) continue;
                        if(count[w] == 0) two_hop_neighbors.push_back(w);
                        count[w]++;
                    }
                }

                Fraction highest_similarity = {0, 1};
                int closest_core = -1;
                for(auto w: two_hop_neighbors){
                    auto sim = Fraction(count[w], g.degree[u] + g.degree[w] - count[w]);
                    auto weight = std::min(sim, similarity_threshold[w]);
                    if((similarity_threshold[w] > similarity_threshold[u]) && (highest_similarity <= weight)){
                        highest_similarity = weight;
                        closest_core = w;
                    }
                    if(visited[w])continue;
                    if(!heap.inside(w)){
                        heap.addToHeap({w, u, std::min(weight, similarity_threshold[u])});
                    }else{
                        heap.increasePriority({w, u, std::min(weight, similarity_threshold[u])});
                    }
                }               
                for(auto w: two_hop_neighbors) count[w] = 0;
                two_hop_neighbors.clear();
                if(closest_core!=-1  && highest_similarity>similarity_threshold[u]) NCI.push_back({closest_core, u,highest_similarity}); //NCI {core, noncore}
            }
            FCF.emplace_back(std::move(currEdges));
            currEdges.clear();
        }
        std::sort(NCI.begin(), NCI.end());
        //need tp sort FCF
    }

    void build_NCI_density(){
        for(auto edges: FCF) FCFTree.push_back(Tree(edges));
        // for(auto edges: FCF) FCFTree1.push_back(Tree1(edges));
        long long allEdges = 0, restEdges = 0;
        std::vector<Edge> edges, compulsoryEdges;
        for(auto vertex: vertices){
            int u = vertex.id;
            for(int i = g.pstart[u]; i<g.pstart[u+1]; i++){
                int v = g.edges[i];
                for(int j = g.pstart[v]; j < g.pstart[v+1]; j++){
                    int w = g.edges[j];
                    if(w==u) continue;
                    if(count[w]==0) two_hop_neighbors.push_back(w);
                    count[w]++;
                }
            }

            
            for(auto w: two_hop_neighbors){
                auto sim = Fraction(count[w], g.degree[u] + g.degree[w] - count[w]);
                auto weight = std::min(sim, similarity_threshold[w]);
                if(weight > similarity_threshold[u]){
                    edges.push_back({w,u,weight});             // NCIEdges {core, noncore}
                }
            }
            sort(edges.begin(), edges.end());
            int treeIdx = rFCF[u];
            auto &tree = FCFTree[treeIdx];
            // auto &tree1 = FCFTree1[treeIdx];
            
            // std::list<EdgeNode> remainEdges;
            for(int i = 0;i<edges.size();i++){
                int v = edges[i].id;
                bool pruned = false;
                // auto it = remainEdges.begin();
                // Fraction currentBestWeight = {0, 1};
                // while(it!=remainEdges.end()){
                //     auto [edge, bestSimilarity] = *it;
                //     if(bestSimilarity >= edges[i].weight){
                //         it = remainEdges.erase(it);
                //         continue;
                //     }
                //     int w = edge.id;
                //     auto minWeight = tree.queryPathMin(tree.getIndex(v), tree.getIndex(w));
                //     // std::cerr << v << " " << w << " " << minWeight << " " << edges[i].weight << std::endl;
                //     currentBestWeight = std::max(minWeight, currentBestWeight);
                //     if(minWeight >= edges[i].weight){
                //         pruned = true;
                //         break;
                //     }
                //     it++;
                // }
                for(int j = 0;j<compulsoryEdges.size();j++){
                    int w = compulsoryEdges[j].id;
                    auto minWeight = tree.queryPathMin(tree.getIndex(v), tree.getIndex(w));
                    // auto minWeight1 = tree1.queryPathMin(tree1.getIndex(v), tree1.getIndex(w));
                    // assert(minWeight1==minWeight);
                    // if(minWeight!=minWeight1){
                        // assert(tree.queryLCA(tree.getIndex(v), tree.getIndex(w))==tree1.lca(tree.getIndex(v), tree.getIndex(w)));
                        // std::cerr << tree.queryLCA(tree.getIndex(v), tree.getIndex(w)) << " " << tree1.lca(tree.getIndex(v), tree.getIndex(w)) << std::endl;
                    //     std::cerr << v << " " << w << " " << minWeight << " " << minWeight1 << std::endl;
                    // }
                    // assert(minWeight==minWeight1);
                    if(minWeight >= edges[i].weight){
                        pruned = true;
                        break;
                    }
                }
                if(!pruned) {
                    compulsoryEdges.push_back(edges[i]);
                    // remainEdges.push_back({edges[i], currentBestWeight});
                }
            }
            allEdges += edges.size(), restEdges += compulsoryEdges.size();
            // std::cerr << u << " " << (double)compulsoryEdges.size() / edges.size() << std::endl;
            for(auto w: two_hop_neighbors) count[w] = 0;
            two_hop_neighbors.clear();
            NCIEdges.insert(NCIEdges.end(),
                std::make_move_iterator(compulsoryEdges.begin()),
                std::make_move_iterator(compulsoryEdges.end()));
            compulsoryEdges.clear();
            edges.clear();
        }
        std::sort(NCIEdges.begin(), NCIEdges.end());
        // std::cout << restEdges << std::endl;
        // std::cout << restEdges << " " << allEdges << std::endl;
        // std::cout << NCI.size() << std::endl;
    }


    void count_NCI_density_full(){
        long long allEdges = 0;
        for(auto vertex: vertices){
            int u = vertex.id;
            for(int i = g.pstart[u]; i<g.pstart[u+1]; i++){
                int v = g.edges[i];
                for(int j = g.pstart[v]; j < g.pstart[v+1]; j++){
                    int w = g.edges[j];
                    if(w==u) continue;
                    if(count[w]==0) two_hop_neighbors.push_back(w);
                    count[w]++;
                }
            }

            
            for(auto w: two_hop_neighbors){
                auto sim = Fraction(count[w], g.degree[u] + g.degree[w] - count[w]);
                auto weight = std::min(sim, similarity_threshold[w]);
                if(weight > similarity_threshold[u]){
                    allEdges++;
                }
            }
            for(auto w: two_hop_neighbors) count[w] = 0;
            two_hop_neighbors.clear();
        }
        std::cerr << "allEdges: " << allEdges << std::endl;
    }

    // std::vector<Cluster> exactClustering(const std::string &directory_name, int mu, const Fraction& eps){
    //     FCFEdges = read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"FCI.bin", eps);
    //     NCI = read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"NCI.bin", eps);
    //     std::vector<Cluster> result;
    //     DisJointSet clusters = DisJointSet(n1);
    //     std::unordered_map<int, int> label2id;
    //     for(auto &edge: FCFEdges){
    //         int u = edge.id, v = edge.id2;
    //         clusters.my_union(u, v);
    //     }

    //     for(auto vertex: vertices){
    //         if(vertex.similarity < eps) break;
    //         int pv = clusters.find_root(vertex.id);
    //         if(!label2id.count(pv)){
    //             label2id[pv] = result.size();
    //             result.push_back(Cluster());
    //         }
    //         result[label2id[pv]].cores.push_back(vertex.id);
    //     }

    //     for(auto edge: NCI){
    //         int v = edge.id, u = edge.id2;  // core, noncore
    //         auto weight = edge.weight; 
    //         if(similarity_threshold[u]<eps){
    //             int pv = clusters.find_root(v);
    //             result[label2id[pv]].nonCores.push_back(u);
    //         }       
    //     }
    //     int count = 0;
    //     for(auto &cluster: result){
    //         // std::cout << cluster.cores.size() << " " << cluster.nonCores.size() << std::endl;
    //         if(cluster.cores.size()==1 && cluster.nonCores.size()==0)count++;
    //     }
    //     std::ofstream log_file("clustering_info.txt", std::ios::app);
    //     std::string msg = "Exact | Dataset " + directory_name + ", mu: " + std::to_string(mu) +
    //                     ", eps: " + std::to_string(eps.num()) + "/" + std::to_string(eps.den()) +
    //                     ", info: " + std::to_string(count) + "/" + std::to_string(result.size());
    //     std::cout << msg << std::endl;
    //     log_file << msg << std::endl;
    //     log_file.close();        
    //     return result;
    // }


    // std::vector<Cluster> densityClustering(const std::string &directory_name, int mu, const Fraction& eps){
    //     FCFEdges = read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"FCI.bin", eps);
    //     NCIEdges = read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"DNCI.bin", eps);
    //     std::vector<Cluster> result;
    //     DisJointSet clusters = DisJointSet(n1);
    //     std::unordered_map<int, int> label2id;
    //     for(auto &edge: FCFEdges){
    //         int u = edge.id, v = edge.id2;
    //         clusters.my_union(u, v);
    //     }

    //     for(auto vertex: vertices){
    //         if(vertex.similarity < eps) break;
    //         int pv = clusters.find_root(vertex.id);
    //         if(!label2id.count(pv)){
    //             label2id[pv] = result.size();
    //             result.push_back(Cluster());
    //         }
    //         result[label2id[pv]].cores.push_back(vertex.id);
    //     }

    //     for(auto edge: NCIEdges){
    //         int v = edge.id, u = edge.id2;  // core, noncore
    //         auto weight = edge.weight; 
    //         if(similarity_threshold[u]<eps){
    //             int pv = clusters.find_root(v);
    //             result[label2id[pv]].nonCores.push_back(u);
    //         }       
    //     }

    //     int count = 0;
    //     for(auto &cluster: result){
            
    //         // std::cout << cluster.cores.size() << " " << cluster.nonCores.size() << std::endl;
    //         if(cluster.cores.size()==1 && cluster.nonCores.size()==0){
    //             count++;
    //         }
    //     }
    //     std::ofstream log_file("clustering_info.txt", std::ios::app);
    //     std::string msg = "Density | Dataset " + directory_name + ", mu: " + std::to_string(mu) +
    //                     ", eps: " + std::to_string(eps.num()) + "/" + std::to_string(eps.den()) +
    //                     ", info: " + std::to_string(count) + "/" + std::to_string(result.size());
    //     std::cout << msg << std::endl;
    //     log_file << msg << std::endl;
    //     log_file.close();        
    //     return result;
    // }


    std::vector<Cluster> exactClustering(const std::string &directory_name, int mu, const Fraction& eps){
        FCFEdges = read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"FCI.bin", eps);
        NCIEdges = read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"NCI.bin", eps);
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

            //     auto elapsed = t.elapsed();
            // std::string msg = "Exact | Dataset: " + directory_name + ", mu: " + std::to_string(mu) +
            //                 ", eps: " + std::to_string(eps.num()) + "/" + std::to_string(eps.den()) +
            //                 ", time(ms): " + Utility::integer_to_string(elapsed);
            // std::cout << msg << std::endl;
            // log_file << msg << std::endl;
            // log_file.close();

        return result;
    }

    std::vector<Cluster> densityClustering(const std::string &directory_name, int mu, const Fraction& eps){
        FCFEdges = std::move(read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"FCI.bin", eps));
        NCIEdges = std::move(read_edges("index/"+directory_name+"_"+std::to_string(mu)+"/"+"DNCI.bin", eps));
        // std::ofstream log_file("clustering_log.txt", std::ios::app);
        // Timer t;
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
            // auto elapsed = t.elapsed();
            // std::string msg = "Density | Dataset: " + directory_name + ", mu: " + std::to_string(mu) +
            //                 ", eps: " + std::to_string(eps.num()) + "/" + std::to_string(eps.den()) +
            //                 ", time(ms): " + Utility::integer_to_string(elapsed);
            // std::cout << msg << std::endl;
            // log_file << msg << std::endl;
            // log_file.close();
        return result;
    }
};