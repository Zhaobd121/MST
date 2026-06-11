#pragma once
#include <vector>
#include "edge.h"
#include <unordered_map>
#include <cassert>
#include <algorithm>

class Tree1{
public:
    int nLabels = 0, n = 0;
    int LOG = 0;
    std::vector<int> h, ne, e, depth;
    std::vector<Fraction> weights;
    std::vector<std::vector<int>> fa;
    std::vector<std::vector<Fraction>> minSimilarity; 
    int idx = 0;  
    std::unordered_map<int,int> labelId;
    std::vector<int> vLabel;

    Tree1(const std::vector<Edge> &edges){
        int m = edges.size();
        if(m>0){
            h.assign(m+1, -1);
            for(auto edge: edges){
                addEdgeWithLabel(edge.id, edge.id2, edge.weight);
                addEdgeWithLabel(edge.id2, edge.id, edge.weight);
            }
            assert(n==m+1);
            buildLCA(0);
        }
    }

    void dfs(int u, int p, int d, Fraction wUp) {
        fa[0][u]     = p;
        minSimilarity[0][u] = wUp;
        depth[u]     = d;

        for (int i = h[u]; i != -1; i = ne[i]) {
            int v = e[i];
            if (v == p) continue;
            dfs(v, u, d + 1, weights[i]);
        }
    }

    void buildLCA(int root) {
        depth.assign(n, 0);
        
        LOG = 1;
        while ((1 << LOG) <= n) ++LOG;
        fa.assign(LOG, std::vector<int>(n, -1));
        minSimilarity.assign(LOG, std::vector<Fraction>(n));

        dfs(root, -1, 0, {2,1});
        

        for (int k = 1; k < LOG; ++k) {
            for (int v = 0; v < n; ++v) {
                int mid = fa[k - 1][v];
                if (mid != -1) {
                    fa[k][v]     = fa[k - 1][mid];
                    minSimilarity[k][v] = std::min(minSimilarity[k - 1][v], minSimilarity[k - 1][mid]);
                }
            }
        }
    }

    int lca(int u, int v) const {
        if (depth[u] < depth[v]) std::swap(u, v);
        int diff = depth[u] - depth[v];
        for (int k = LOG - 1; k >= 0; --k)
            if (diff & (1 << k))
                u = fa[k][u];

        if (u == v) return u;

        for (int k = LOG - 1; k >= 0; --k)
            if (fa[k][u] != fa[k][v]) {
                u = fa[k][u];
                v = fa[k][v];
            }
        return fa[0][u];
    }

    Fraction queryPathMin(int u, int v) const {
        Fraction ans = {2, 1};
        if (depth[u] < depth[v]) std::swap(u, v);

        int diff = depth[u] - depth[v];
        for (int k = LOG - 1; k >= 0; --k)
            if (diff & (1 << k)) {
                ans = std::min(ans, minSimilarity[k][u]);
                u   = fa[k][u];
            }

        if (u == v) return ans;

        for (int k = LOG - 1; k >= 0; --k)
            if (fa[k][u] != fa[k][v]) {
                ans = std::min({ans, minSimilarity[k][u], minSimilarity[k][v]});
                u = fa[k][u];
                v = fa[k][v];
            }
        ans = std::min({ans, minSimilarity[0][u], minSimilarity[0][v]});
        return ans;
    }

    int getIndex(int label) {
        if (labelId.count(label)) return labelId[label];
        vLabel.push_back(label);
        return labelId[label] = n++;
    }

    void addEdgeWithLabel(int ulabel, int vlabel, Fraction weight) {
        nLabels = std::max(nLabels, ulabel + 1);
        nLabels = std::max(nLabels, vlabel + 1);
        addEdge(getIndex(ulabel), getIndex(vlabel), weight);
    }

    void addEdge(int u,int v,Fraction wt) {
        e.push_back(v);
        weights.push_back(wt);
        ne.push_back(h[u]);
        h[u] = e.size() - 1;
    }
};