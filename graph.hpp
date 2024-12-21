//Author : Ruben Rehal
//Date : November 30, 2024
//Advanced Graphing

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>

// Exception class definition
class illegal_exception {
public:
    const char* what() const { return "illegal argument"; }
};

// Edge Class
class Edge {
public:
    Edge(double w, const std::string& l);
    
    double getWeight() const;
    void setWeight(double w);
    
    std::string getLabel() const;
    void setLabel(const std::string& l);
    
private:
    double weight;
    std::string label;
};

// Node Class
class Node {
public:
    Node(const std::string& l, const std::string& ID, const std::string& T);
    
    std::string getID() const;
    std::string getLabel() const;
    void setLabel(const std::string& l);
    
    std::string getType() const;
    void setType(const std::string& t);
    
private:
    std::string label; // Name
    std::string id;
    std::string type;
};

// Adjacency List Node Class
class AdjacencyNode {
public:
    AdjacencyNode(int neighbor, Edge* e, AdjacencyNode* n = nullptr);
    
    int getNeighborIndex() const;
    Edge* getEdge() const;
    AdjacencyNode* getNext() const;
    void setNext(AdjacencyNode* n);
    
private:
    int neighborIndex;
    Edge* edge;
    AdjacencyNode* next;
};

// Graph Class
class Graph {
public:
    Graph();
    ~Graph();
    
    void load(const std::string& filename, const std::string& type);
    void addNode(const std::string& id, const std::string& name, const std::string& type);
    bool addRelationship(const std::string& sourceID, const std::string& label, const std::string& destID, double weight);
    void deleteNode(const std::string& nodeID);
    void printAdjacencyList(const std::string& vertexID);
    void path(const std::string& ID_1, const std::string& ID_2);
    void findMaxWeightPath();
    void findAll(const std::string& fieldType, const std::string& fieldString);
    
private:
    void resize();
    void addEdge(int node1Index, int node2Index, double weight, const std::string& edgeLabel);
    bool hasEdges();
    
    Node** nodes;                 // Dynamic array of node pointers
    AdjacencyNode** adjLists;     // Dynamic array of adjacency list heads
    int nodeCount;
    int capacity;
};

#endif // GRAPH_HPP
