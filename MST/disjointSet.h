#pragma once
#include <vector>

class DisJointSet{
public:
    std::vector<int> parent, rank;

    DisJointSet(int n): parent(std::vector<int>(n)), rank(std::vector<int>(n, 0)){
        for(int i = 0;i<n;i++) parent[i] = i;
    }

    int find_root(int u) {
        int x = u;
        while(parent[x] != x) x = parent[x];

        while(parent[u] != x) {
            int tmp = parent[u];
            parent[u] = x;
            u = tmp;
        }
        return x;
    }

    void my_union(int u, int v) {
        int pu = find_root(u);
        int pv = find_root(v);

        if(pu == pv) return ;

        if(rank[pu] < rank[pv]) parent[pu] = pv;
        else if(rank[pu] > rank[pv]) parent[pv] = pu;
        else {
                parent[pv] = pu;
                ++ rank[pu];
        }
    }
};