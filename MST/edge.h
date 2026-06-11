#pragma once
#include "Fraction.h"

class Edge{
public:
    int id, id2;
    Fraction weight;

    Edge() = default;

    Edge(int u, int v, Fraction w)                       
        : id(u), id2(v), weight(std::move(w)) {}

    bool operator<(const Edge& other) const {
        if(weight != other.weight) return weight > other.weight;  
        return id < other.id;
    }
};