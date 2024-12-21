//Author : Ruben Rehal
//Date : November 30, 2024
//Advanced Graphing

#include "graph.hpp"

// Edge Class Implementation
Edge::Edge(double w, const std::string& l) : weight(w), label(l) {}

double Edge::getWeight() const {
    return weight;
}

void Edge::setWeight(double w) {
    weight = w;
}

std::string Edge::getLabel() const {
    return label;
}

void Edge::setLabel(const std::string& l) {
    label = l;
}

// Node Class Implementation
Node::Node(const std::string& l, const std::string& ID, const std::string& T) : label(l), id(ID), type(T) {}

std::string Node::getID() const {
    return id;
}

std::string Node::getLabel() const {
    return label;
}

void Node::setLabel(const std::string& l) {
    label = l;
}

std::string Node::getType() const {
    return type;
}

void Node::setType(const std::string& t) {
    type = t;
}

// AdjacencyNode Class Implementation
AdjacencyNode::AdjacencyNode(int neighbor, Edge* e, AdjacencyNode* n) : neighborIndex(neighbor), edge(e), next(n) {}

int AdjacencyNode::getNeighborIndex() const {
    return neighborIndex;
}

Edge* AdjacencyNode::getEdge() const {
    return edge;
}

AdjacencyNode* AdjacencyNode::getNext() const {
    return next;
}

void AdjacencyNode::setNext(AdjacencyNode* n) {
    next = n;
}

// Graph Class Implementation

// Constructor
Graph::Graph() : nodeCount(0), capacity(2) {
    nodes = new Node*[capacity];
    adjLists = new AdjacencyNode*[capacity];
    for (int i = 0; i < capacity; i++) {
        adjLists[i] = nullptr;
    }
}

// Destructor
Graph::~Graph() {
    // Delete all adjacency lists
    for (int i = 0; i < nodeCount; i++) {
        AdjacencyNode* current = adjLists[i];
        while (current != nullptr) {
            AdjacencyNode* temp = current;
            current = current->getNext();
            // Delete the Edge only once, as edges are stored in both adjacency lists
            // To prevent double deletion, only delete the Edge when neighborIndex > i
            if (temp->getNeighborIndex() > i && temp->getEdge() != nullptr) {
                delete temp->getEdge();
            }
            delete temp;
        }
    }
    // Delete adjacency list array
    delete[] adjLists;
    
    // Delete all nodes
    for (int i = 0; i < nodeCount; i++) {
        delete nodes[i];
    }
    // Delete nodes array
    delete[] nodes;
}

// Resize function to double the capacity
void Graph::resize() {
    int newCapacity = capacity * 2;
    Node** newNodes = new Node*[newCapacity];
    AdjacencyNode** newAdjLists = new AdjacencyNode*[newCapacity];
    
    // Copy existing nodes and adjLists
    for (int i = 0; i < nodeCount; i++) {
        newNodes[i] = nodes[i];
        newAdjLists[i] = adjLists[i];
    }
    
    // Initialize new adjLists
    for (int i = nodeCount; i < newCapacity; i++) {
        newAdjLists[i] = nullptr;
    }
    
    // Delete old arrays
    delete[] nodes;
    delete[] adjLists;
    
    // Update pointers and capacity
    nodes = newNodes;
    adjLists = newAdjLists;
    capacity = newCapacity;
}

// Load function with type
void Graph::load(const std::string& filename, const std::string& type) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "failure" << std::endl;
        return;
    }

    std::string line;
    if (type == "entities") {
        // Load entities: id, name, type
        while (std::getline(file, line)) {
            if (!line.empty()) {
                std::istringstream iss(line);
                std::string id, name, nodeType;
                iss >> id >> name >> nodeType;
                addNode(id, name, nodeType);
            }
        }
    }
    else if (type == "relationships") {
        // Load relationships: sourceID, label, destID, weight
        while (std::getline(file, line)) {
            if (!line.empty()) {
                std::istringstream iss(line);
                std::string sourceID, label, destID;
                double weight;
                iss >> sourceID >> label >> destID >> weight;
                addRelationship(sourceID, label, destID, weight);
            }
        }
    }
    else {
        // Invalid type (should not occur as per project constraints)
        std::cout << "failure" << std::endl;
        file.close();
        return;
    }

    file.close();
    std::cout << "success" << std::endl;
}

// Add a node (entity)
void Graph::addNode(const std::string& id, const std::string& name, const std::string& type) {
    // Check if the node already exists; if so, update its name and type
    for (int i = 0; i < nodeCount; i++) {
        if (nodes[i]->getID() == id) {
            nodes[i]->setLabel(name);
            nodes[i]->setType(type);
            return;
        }
    }

    // If node does not exist, add it
    if (nodeCount >= capacity) {
        resize();
    }

    nodes[nodeCount] = new Node(name, id, type);
    adjLists[nodeCount] = nullptr; // Initialize adjacency list
    nodeCount++;
}

// Add or update a relationship (edge)
bool Graph::addRelationship(const std::string& sourceID, const std::string& label, const std::string& destID, double weight) {
    // If Negative or Zero Weight
    if ( weight <= 0 ) {
        return false;
    }

    // Find the indices of sourceID and destID
    int sourceIndex = -1, destIndex = -1;
    for (int i = 0; i < nodeCount; i++) {
        if (nodes[i]->getID() == sourceID) sourceIndex = i;
        if (nodes[i]->getID() == destID) destIndex = i;
    }

    // If either node doesn't exist, ignore this edge as per project specifications
    if (sourceIndex == -1 || destIndex == -1) {
        return false;
    }

    // Check if an edge already exists between source and dest
    AdjacencyNode* current = adjLists[sourceIndex];
    while (current != nullptr) {
        if (current->getNeighborIndex() == destIndex) {
            // Edge exists, update it
            current->getEdge()->setLabel(label);
            current->getEdge()->setWeight(weight);
            // Also update the reverse edge
            AdjacencyNode* reverse = adjLists[destIndex];
            while (reverse != nullptr) {
                if (reverse->getNeighborIndex() == sourceIndex) {
                    reverse->getEdge()->setLabel(label);
                    reverse->getEdge()->setWeight(weight);
                    break;
                }
                reverse = reverse->getNext();
            }
            return true;
        }
        current = current->getNext();
    }

    // If no existing edge, add a new edge
    addEdge(sourceIndex, destIndex, weight, label);
    return true;
}

// Add an edge to the graph (internal use)
void Graph::addEdge(int node1Index, int node2Index, double weight, const std::string& edgeLabel) {
    // Create a new edge
    Edge* newEdge = new Edge(weight, edgeLabel);
    
    // Add to node1's adjacency list
    AdjacencyNode* newAdj1 = new AdjacencyNode(node2Index, newEdge, adjLists[node1Index]);
    adjLists[node1Index] = newAdj1;
    
    // Add to node2's adjacency list
    AdjacencyNode* newAdj2 = new AdjacencyNode(node1Index, newEdge, adjLists[node2Index]);
    adjLists[node2Index] = newAdj2;
}

// Delete a node by ID
void Graph::deleteNode(const std::string& nodeID) {
    try {
        // Validate nodeID
        for (char c : nodeID) {
            if (!std::isalnum(static_cast<unsigned char>(c))) {
                throw illegal_exception();
            }
        }

        // Find the node index by its ID
        int nodeIndex = -1;
        for (int i = 0; i < nodeCount; i++) {
            if (nodes[i]->getID() == nodeID) {
                nodeIndex = i;
                break;
            }
        }

        // If node is not found, print "failure" and return
        if (nodeIndex == -1) {
            std::cout << "failure" << std::endl;
            return;
        }

        // Remove all edges from this node
        AdjacencyNode* current = adjLists[nodeIndex];
        while (current != nullptr) {
            int neighbor = current->getNeighborIndex();
            // Remove the edge from neighbor's adjacency list
            AdjacencyNode* prevAdj = nullptr;
            AdjacencyNode* adjCurrent = adjLists[neighbor];
            while (adjCurrent != nullptr) {
                if (adjCurrent->getNeighborIndex() == nodeIndex) {
                    // Found the adjacency node to remove
                    if (prevAdj == nullptr) {
                        adjLists[neighbor] = adjCurrent->getNext();
                    }
                    else {
                        prevAdj->setNext(adjCurrent->getNext());
                    }
                    // If neighbor > nodeIndex, delete the Edge
                    if (neighbor > nodeIndex && adjCurrent->getEdge() != nullptr) {
                        delete adjCurrent->getEdge();
                    }
                    delete adjCurrent;
                    break;
                }
                prevAdj = adjCurrent;
                adjCurrent = adjCurrent->getNext();
            }
            // Move to next adjacency node
            AdjacencyNode* temp = current;
            current = current->getNext();
            // Delete the adjacency node from this node's list
            delete temp;
        }
        adjLists[nodeIndex] = nullptr;

        // Delete the node itself
        delete nodes[nodeIndex];

        // Shift nodes and adjacency lists to fill the gap
        for (int i = nodeIndex; i < nodeCount - 1; i++) {
            nodes[i] = nodes[i + 1];
            adjLists[i] = adjLists[i + 1];
            // Update neighbor indices in adjacency lists
            AdjacencyNode* adj = adjLists[i];
            while (adj != nullptr) {
                if (adj->getNeighborIndex() > nodeIndex) {
                    // Decrement neighbor index since nodes have shifted
                    adj->setNext(adj->getNext());
                }
                adj = adj->getNext();
            }
        }

        // Set the last node and adjacency list to nullptr
        nodes[nodeCount - 1] = nullptr;
        adjLists[nodeCount - 1] = nullptr;

        // Decrease the node count
        nodeCount--;

        std::cout << "success" << std::endl;
    }
    catch (const illegal_exception& e) {
        std::cout << e.what() << std::endl;
    }
}

// Print adjacency list of a node
void Graph::printAdjacencyList(const std::string& vertexID) {
    try {
        // Validate vertexID
        for (char c : vertexID) {
            if (!std::isalnum(static_cast<unsigned char>(c))) {
                throw illegal_exception();
            }
        }

        // Find the index of the vertex with the given ID
        int vertexIndex = -1;
        for (int i = 0; i < nodeCount; i++) {
            if (nodes[i]->getID() == vertexID) {
                vertexIndex = i;
                break;
            }
        }

        // If the vertex is not found, print "failure" and return
        if (vertexIndex == -1) {
            std::cout << "failure" << std::endl;
            return;
        }

        // Collect all adjacent vertex IDs
        std::string adjacencyOutput = "";
        AdjacencyNode* current = adjLists[vertexIndex];
        while (current != nullptr) {
            if (!adjacencyOutput.empty()) {
                adjacencyOutput += " ";
            }
            adjacencyOutput += nodes[current->getNeighborIndex()]->getID();
            current = current->getNext();
        }

        // Print the adjacency list or a blank line
        std::cout << adjacencyOutput << std::endl;
    }
    catch (const illegal_exception& e) {
        std::cout << e.what() << std::endl;
    }
}

// Path function using Dijkstra's algorithm variant for maximum weight paths
void Graph::path(const std::string& ID_1, const std::string& ID_2) {
    try {
        // Validate IDs
        for (char c : ID_1) {
            if (!std::isalnum(static_cast<unsigned char>(c))) {
                throw illegal_exception();
            }
        }
        for (char c : ID_2) {
            if (!std::isalnum(static_cast<unsigned char>(c))) {
                throw illegal_exception();
            }
        }

        // Check if the graph is empty
        if (nodeCount == 0) {
            std::cout << "failure" << std::endl;
            return;
        }

        // Check if both vertices exist in the graph
        int startIdx = -1, endIdx = -1;
        for (int i = 0; i < nodeCount; i++) {
            if (nodes[i]->getID() == ID_1) {
                startIdx = i;
            }
            if (nodes[i]->getID() == ID_2) {
                endIdx = i;
            }
        }

        // If either vertex is not found, print "failure"
        if (startIdx == -1 || endIdx == -1) {
            std::cout << "failure" << std::endl;
            return;
        }

        // Initialize distance array and predecessor array for the path reconstruction
        double* dist = new double[nodeCount];
        int* prev = new int[nodeCount];
        bool* visited = new bool[nodeCount];
        
        for (int i = 0; i < nodeCount; i++) {
            dist[i] = -INFINITY; // Start with the most negative number (since we're looking for the highest weight)
            prev[i] = -1; // No predecessors initially
            visited[i] = false; // No nodes visited yet
        }

        dist[startIdx] = 0; // Starting vertex has a distance of 0

        // Simple greedy approach to find the highest weight path (modified Dijkstra's)
        for (int i = 0; i < nodeCount; i++) {
            // Find the unvisited node with the highest distance
            int current = -1;
            double maxDist = -INFINITY;

            for (int j = 0; j < nodeCount; j++) {
                if (!visited[j] && dist[j] > maxDist) {
                    maxDist = dist[j];
                    current = j;
                }
            }

            if (current == -1) break; // If no more reachable node, break out

            visited[current] = true;

            // Traverse adjacency list of current
            AdjacencyNode* neighbor = adjLists[current];
            while (neighbor != nullptr) {
                int neighborIdx = neighbor->getNeighborIndex();
                double edgeWeight = neighbor->getEdge()->getWeight();

                if (!visited[neighborIdx]) {
                    double newDist = dist[current] + edgeWeight;
                    if (newDist > dist[neighborIdx]) {
                        dist[neighborIdx] = newDist;
                        prev[neighborIdx] = current;
                    }
                }
                neighbor = neighbor->getNext();
            }
        }

        // If there's no path to the end node
        if (dist[endIdx] == -INFINITY) {
            std::cout << "failure" << std::endl;
            delete[] dist;
            delete[] prev;
            delete[] visited;
            return;
        }

        // Reconstruct the path
        std::string pathString = "";
        int currentNode = endIdx;
        double totalWeight = dist[endIdx];
        
        while (currentNode != -1) {
            pathString = nodes[currentNode]->getID() + " " + pathString;
            currentNode = prev[currentNode];
        }

        // Remove the trailing space
        if (!pathString.empty() && pathString.back() == ' ') {
            pathString.pop_back();
        }

        // Manually format the totalWeight to one decimal place
        double roundedWeight = floor(totalWeight * 10 + 0.5) / 10.0;
        int integerPart = static_cast<int>(roundedWeight);
        int decimalPart = static_cast<int>((roundedWeight - integerPart) * 10 + 0.5);

        // Ensure that decimalPart is between 0 and 9
        if (decimalPart == 10) {
            integerPart += 1;
            decimalPart = 0;
        }

        // Print the path and the weight with one decimal place
        std::cout << pathString << " " << integerPart << "." << decimalPart << std::endl;

        // Clean up
        delete[] dist;
        delete[] prev;
        delete[] visited;
    }
    catch (const illegal_exception& e) {
        std::cout << e.what() << std::endl;
    }
}

// Find the two nodes with the highest weight path between them
void Graph::findMaxWeightPath() {
    // Check if the graph is empty or has no edges
    if (nodeCount == 0 || !hasEdges()) {
        std::cout << "failure" << std::endl;
        return;
    }

    double maxWeight = -INFINITY;  // Start with the most negative number
    std::string ID_1 = "", ID_2 = "";
    
    // Loop through all unique pairs of nodes to find the maximum weight path
    for (int i = 0; i < nodeCount; i++) {
        for (int j = i + 1; j < nodeCount; j++) {
            // Perform a path search between nodes[i] and nodes[j]
            // Initialize distance array and predecessor array
            double* dist = new double[nodeCount];
            int* prev = new int[nodeCount];
            bool* visited = new bool[nodeCount];

            for (int k = 0; k < nodeCount; k++) {
                dist[k] = -INFINITY; // Set distance to negative infinity
                prev[k] = -1; // No predecessors initially
                visited[k] = false; // No nodes visited yet
            }

            dist[i] = 0; // Start node has distance 0

            // Greedy approach to find the highest weight path (modified Dijkstra's)
            for (int k = 0; k < nodeCount; k++) {
                // Find the unvisited node with the highest distance
                int currentNode = -1;
                double maxDist = -INFINITY;

                for (int l = 0; l < nodeCount; l++) {
                    if (!visited[l] && dist[l] > maxDist) {
                        maxDist = dist[l];
                        currentNode = l;
                    }
                }

                if (currentNode == -1) break; // No more reachable nodes

                visited[currentNode] = true;

                // Traverse adjacency list of current
                AdjacencyNode* neighbor = adjLists[currentNode];
                while (neighbor != nullptr) {
                    int neighborIdx = neighbor->getNeighborIndex();
                    double edgeWeight = neighbor->getEdge()->getWeight();

                    if (!visited[neighborIdx]) {
                        double newDist = dist[currentNode] + edgeWeight;
                        if (newDist > dist[neighborIdx]) {
                            dist[neighborIdx] = newDist;
                            prev[neighborIdx] = currentNode;
                        }
                    }
                    neighbor = neighbor->getNext();
                }
            }

            // Check if the path from node i to node j has a weight greater than maxWeight
            if (dist[j] > maxWeight) {
                maxWeight = dist[j];
                ID_1 = nodes[i]->getID();
                ID_2 = nodes[j]->getID();
            }

            // Clean up
            delete[] dist;
            delete[] prev;
            delete[] visited;
        }
    }

    // If no valid path is found, print "failure"
    if (maxWeight == -INFINITY) {
        std::cout << "failure" << std::endl;
    }
    else {
        // Manually format the maxWeight to one decimal place
        double roundedWeight = floor(maxWeight * 10 + 0.5) / 10.0;
        int integerPart = static_cast<int>(roundedWeight);
        int decimalPart = static_cast<int>((roundedWeight - integerPart) * 10 + 0.5);

        // Ensure that decimalPart is between 0 and 9
        if (decimalPart == 10) {
            integerPart += 1;
            decimalPart = 0;
        }

        std::cout << ID_1 << " " << ID_2 << " " << integerPart << "." << decimalPart << std::endl;
    }

}

// Check if the graph has any edges
bool Graph::hasEdges() {
    for (int i = 0; i < nodeCount; i++) {
        if (adjLists[i] != nullptr) {
            return true;
        }
    }
    return false;
}

// Find all nodes matching a field
void Graph::findAll(const std::string& fieldType, const std::string& fieldString) {
    bool found = false;
    std::string findAllOutput = "";
    
    // Loop through all nodes
    for (int i = 0; i < nodeCount; i++) {
        // If Field_type is "name", check the label
        if (fieldType == "name" && nodes[i]->getLabel() == fieldString) {
            if (found) {
                findAllOutput += " ";
            }
            findAllOutput += nodes[i]->getID();
            found = true;
        }
        // If Field_type is "type", check the type
        else if (fieldType == "type" && nodes[i]->getType() == fieldString) {
            if (found) {
                findAllOutput += " ";
            }
            findAllOutput += nodes[i]->getID();
            found = true;
        }
    }

    if (!found) {
        std::cout << "failure" << std::endl;
    }
    else {
        std::cout << findAllOutput << std::endl;
    }
}
