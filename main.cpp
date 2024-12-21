// Author : Ruben Rehal
// Date : November 30, 2024
// Advanced Graphing

#include "graph.hpp"

int main() {
    Graph graph;
    std::string command;
    
    while (std::cin >> command) {
        if (command == "LOAD") {
            std::string filename, type;
            std::cin >> filename >> type;
            graph.load(filename, type);
        }
        else if (command == "ENTITY") {
            std::string id, name, type;
            std::cin >> id >> name >> type;
            graph.addNode(id, name, type);
            std::cout << "success" << std::endl;
        }
        else if (command == "RELATIONSHIP") {
            std::string sourceID, label, destID;
            double weight;
            std::cin >> sourceID >> label >> destID >> weight;
            bool result = graph.addRelationship(sourceID, label, destID, weight);
            if (result) {
                std::cout << "success" << std::endl;
            }
            else {
                std::cout << "failure" << std::endl;
            }
        }
        else if (command == "PRINT") {
            std::string vertexID;
            std::cin >> vertexID;
            graph.printAdjacencyList(vertexID);
        }
        else if (command == "DELETE") {
            std::string nodeID;
            std::cin >> nodeID;
            graph.deleteNode(nodeID);
        }
        else if (command == "PATH") {
            std::string ID_1, ID_2;
            std::cin >> ID_1 >> ID_2;
            graph.path(ID_1, ID_2);
        }
        else if (command == "HIGHEST") {
            graph.findMaxWeightPath();
        }
        else if (command == "FINDALL") {
            std::string fieldType, fieldString;
            std::cin >> fieldType >> fieldString;
            graph.findAll(fieldType, fieldString);
        }
        else if (command == "EXIT") {
            break;
        }
        // Ignore any other commands as per project specifications
    }
    
    return 0;
}
