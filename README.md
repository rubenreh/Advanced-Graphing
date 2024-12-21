# Advanced Graph Processing System

## Overview

This repository contains a sophisticated implementation of an **Advanced Graph Processing System**, designed to tackle the most demanding computational graph challenges. Leveraging cutting-edge algorithms and efficient memory management, this project redefines graph processing capabilities.

---

## Key Features

### 1. **Dynamic Graph Representation**
- **Scalable Design**: Dynamic resizing of nodes and adjacency lists ensures efficient use of memory while accommodating extensive graph expansions.
- **Bi-directional Relationships**: Each edge is seamlessly integrated into the adjacency lists of both source and destination nodes for bidirectional traversal.

### 2. **Algorithmic Innovations**
- **Modified Dijkstra's Algorithm**: An enhanced version that prioritizes maximum weight paths, ensuring optimal traversal strategies.
- **Custom Pathfinding**: Implementations for finding the maximum weight path between two nodes or globally across all nodes.
- **Entity and Relationship Loading**: Batch processing of nodes and edges from files, offering rapid initialization.

### 3. **Robust Memory Management**
- **Custom Destructor**: Handles the dynamic cleanup of nodes, edges, and adjacency lists to prevent memory leaks.
- **Efficient Resizing**: Employs exponential growth strategies to minimize resizing overhead.

### 4. **Error Handling**
- **Validation Mechanisms**: Stringent checks for illegal inputs, invalid node IDs, and graph inconsistencies.
- **Graceful Failures**: Clear and concise failure responses for invalid operations.

### 5. **Extensive Node and Edge Operations**
- **Dynamic Node Addition**: Handles duplication gracefully with updates to labels and types.
- **Relationship Updates**: Existing edges are updated dynamically without redundant allocations.
- **Field-Based Queries**: Search nodes by attributes (`name` or `type`) with highly optimized filters.

---

## Technical Implementation

### Core Classes

#### `Node`
Encapsulates graph vertices with attributes:
- `label`: User-defined descriptor.
- `id`: Unique identifier.
- `type`: Node classification.

#### `Edge`
Manages connections between nodes, storing:
- `weight`: Represents the cost or significance of the connection.
- `label`: Descriptive metadata.

#### `AdjacencyNode`
Forms the backbone of adjacency lists:
- Links neighboring nodes via edge references.
- Allows seamless graph traversal.

#### `Graph`
Orchestrates the system with:
- **Dynamic Memory Handling**: Implements a dual-pointer system for nodes and adjacency lists.
- **Algorithms**: Houses pathfinding, maximum path weight discovery, and global edge analysis.

### Algorithms
#### Modified Dijkstra's Algorithm
This variant is tailored for applications prioritizing **maximum weight paths**. It balances speed and computational accuracy by:
- Using a priority queue structure (simulated in this implementation).
- Computing cumulative maximum weights instead of minimum distances.

#### Graph Expansion and Optimization
- Dynamic allocation with doubling strategies ensures **O(log N)** memory growth.
- Efficient node and edge lookup operations are enabled by custom index mapping.

---

## Usage

### Compilation
```bash
g++ -std=c++17 -o graph_processor main.cpp p4.cpp
```

### Running the System
```bash
./graph_processor <command_file>
```

### Supported Commands
1. **Load Graph**
   ```text
   load <filename> <type>  # type: entities or relationships
   ```
2. **Add Node**
   ```text
   addNode <id> <name> <type>
   ```
3. **Add Relationship**
   ```text
   addRelationship <sourceID> <label> <destID> <weight>
   ```
4. **Query Path**
   ```text
   path <ID_1> <ID_2>
   ```
5. **Find Maximum Weight Path**
   ```text
   findMaxWeightPath
   ```

---

## Future Enhancements
- **Graph Visualization**: Integration with tools like Graphviz for intuitive visual representations.
- **Asynchronous Processing**: Employing parallel algorithms for faster execution on large graphs.
- **Persistent Storage**: Incorporating a database backend for real-time updates and querying.

---

## Citations
Special thanks to the contributors of graph theory literature, whose insights significantly influenced the **pathfinding and optimization algorithms** in this project.
