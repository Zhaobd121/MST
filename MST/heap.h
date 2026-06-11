#pragma once
#include <vector>
#include <cassert>
#include "edge.h"

class Heap{
public:
    std::vector<Edge> heap;
    std::vector<int> pos;
    int heap_n;

    Heap(int n): heap(std::vector<Edge>(n)), pos(std::vector<int>(n, -1)) , heap_n(0){}

    bool inside(int idx){return pos[idx]!=-1;}

    int size(){return heap_n;}

    void increasePriority(const Edge &node){
        assert(inside(node.id));
        auto tmp = heap[pos[node.id]];
        if(node<tmp){
            assert(node.weight>=tmp.weight);
            heap[pos[node.id]] = node;
            heapBottomUp(pos[node.id]);
        }
    }

    void addToHeap(const Edge &node){
        assert(pos[node.id] == -1);
        heap[heap_n] = node;
        pos[node.id] = heap_n;
        ++ heap_n;
        heapBottomUp(heap_n-1);
    }

    void heapBottomUp(int idx) {
        auto tmp = heap[idx];
        while(idx > 0) {
            int i = (idx-1)/2;
            if(tmp < heap[i]) {
                heap[idx] = heap[i];
                pos[heap[idx].id] = idx;
                idx = i;
            }
            else break;
        }
        heap[idx] = tmp;
        pos[tmp.id] = idx;        
    }

    void heapTopDown(int idx) {
        auto tmp = heap[idx];
        while(2*idx+1 < heap_n) {
            int i = 2*idx+1;
            if(i+1 < heap_n&&heap[i+1] < heap[i]) ++ i;
            if(heap[i] < tmp) {
                heap[idx] = heap[i];
                pos[heap[idx].id] = idx;
                idx = i;
            }
            else break;
        }
        heap[idx] = tmp;
        pos[tmp.id] = idx;
    }

    Edge pop(){
        auto res = heap[0];
        pos[res.id] = -1;
        heap_n--;
        if(heap_n){
            heap[0] = heap[heap_n];
            pos[heap[0].id] = 0;
            heapTopDown(0);
        }
        return res;
    }

};