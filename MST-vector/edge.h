#pragma once

class Edge{
public:
    int id, id2;
    double weight;

    Edge() = default;

    Edge(int u, int v, double w)                       
        : id(u), id2(v), weight(w) {}

    bool operator<(const Edge& other) const {
        if(weight != other.weight) return weight < other.weight;  
        return id < other.id;
    }
};