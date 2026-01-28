#pragma once
#include <vector>
#include "edge.h"
#include <unordered_map>
#include <cassert>
#include <algorithm>

class Tree{
public:
    int nLabels = 0, n = 0;
    int LOG = 0;
    const double INF = 9999999;
    std::vector<int> h, ne, e;
    std::vector<double> weights;
    int idx = 0;  
    std::unordered_map<int,int> labelId;
    std::vector<int> vLabel;
    

    Tree(std::vector<Edge> &edges){
        int m = edges.size();
        if(m>0){
            h.assign(2*m+1, -1);
            DisJointSet disJointSet = DisJointSet(m+1);
            std::vector<int> rep(m+1);
            for(int i = 0;i<m+1;i++)rep[i] = i;
            for(auto edge: edges){
                getIndex(edge.id);
                getIndex(edge.id2);
            }
            std::sort(edges.begin(), edges.end());
            for(auto edge: edges){
                int u = getIndex(edge.id), v = getIndex(edge.id2);
                int pu = disJointSet.find_root(u), pv = disJointSet.find_root(v);
                if(pu==pv)continue;
                weights.push_back(edge.weight);
                addEdge(rep[pu], n), addEdge(n, rep[pu]);
                addEdge(rep[pv], n), addEdge(n, rep[pv]);
                disJointSet.my_union(u, v);
                rep[disJointSet.find_root(u)] = n++;
            }
            assert(weights.size()==n);
            assert(weights.size()==2*m+1);
            assert(e.size()==4*m);
            build(n-1);
        }
    }

    std::vector<int> depth_, firstOcc_, euler_, depthEuler_;
    std::vector<int> lg_;
    std::vector<std::vector<int>> stDepth_;  

    void dfs(int u,int p){
        firstOcc_[u] = (int)euler_.size();
        euler_.push_back(u);
        depthEuler_.push_back(depth_[u]);

        for(int i=h[u]; i!=-1; i=ne[i]){
            int v = e[i];
            if(v==p) continue;
            depth_[v] = depth_[u] + 1;
            dfs(v,u);
            euler_.push_back(u);
            depthEuler_.push_back(depth_[u]);
        }
    }  

    void build(int root){
        depth_.assign(n,0);
        firstOcc_.assign(n,-1);
        euler_.reserve(2*n);
        depthEuler_.reserve(2*n);
        dfs(root,-1);

        int m = euler_.size();
        lg_.assign(m+1,0);
        for(int i=2;i<=m;++i) lg_[i] = lg_[i>>1]+1;
        int K = lg_[m]+1;

        stDepth_.assign(K, std::vector<int>(m));
        for(int i=0;i<m;++i) stDepth_[0][i] = i;
        for(int k=1, len=1; k<K; ++k, len<<=1)
            for(int i=0; i+(len<<1)<=m; ++i){
                int a = stDepth_[k-1][i];
                int b = stDepth_[k-1][i+len];
                stDepth_[k][i] =
                    depthEuler_[a] < depthEuler_[b] ? a : b;
            }
    }

    int lca(int u,int v) const{
        int l = firstOcc_[u], r = firstOcc_[v];
        if(l>r) std::swap(l,r);
        int k = lg_[r-l+1];
        int a = stDepth_[k][l];
        int b = stDepth_[k][r-(1<<k)+1];
        return depthEuler_[a] < depthEuler_[b] ? euler_[a] : euler_[b];
    }

    double queryPathMin(int u, int v) const{
        return weights[lca(u, v)];
    }

    int getIndex(int label) {
        if (labelId.count(label)) return labelId[label];
        vLabel.push_back(label);
        weights.push_back(0);
        return labelId[label] = n++;
    }

    void addEdgeWithLabel(int ulabel, int vlabel) {
        nLabels = std::max(nLabels, ulabel + 1);
        nLabels = std::max(nLabels, vlabel + 1);
        addEdge(getIndex(ulabel), getIndex(vlabel));
    }

    void addEdge(int u,int v) {
        e.push_back(v);
        ne.push_back(h[u]);
        h[u] = e.size() - 1;
    }
};