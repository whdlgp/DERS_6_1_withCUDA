#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>

namespace PushRelabel
{
    
struct Edge
{
	// To store current flow and capacity of edge 
	int flow, capacity;

	// An edge u--->v has start vertex as u and end 
	// vertex as v. 
	int u, v;

	Edge(int flow, int capacity, int u, int v)
	{
		this->flow = flow;
		this->capacity = capacity;
		this->u = u;
		this->v = v;
	}
};

// Represent a Vertex 
struct Vertex
{
	int h, e_flow;

	Vertex(int h, int e_flow)
	{
		this->h = h;
		this->e_flow = e_flow;
	}
};

// To represent a flow network 
class Graph
{
	int V;    // No. of vertices 
	std::vector<Vertex> ver;
	std::vector<Edge> edge;
	std::vector<std::vector<int>> edge_map; // Hashmap of edge, find v with key u
	std::vector<int> e_flow_map;

	std::vector<bool> visited;

	int s;
	int t;

	// Function to push excess flow from u 
	bool push(int u);

	// Function to relabel a vertex u 
	void relabel(int u);

	// This function is called to initialize 
	// preflow 
	void preflow(int s);

	// Function to reverse edge 
	void updateReverseEdgeFlow(int i, int flow);

public:
	Graph(int V);  // Constructor
    Graph();

				   // function to add an edge to graph 
	void addEdge(int u, int v, int w);
                   // function to add a node to graph
    int addNode();

	               // returns maximum flow from s to t 
	int getMaxFlow(int s, int t);
                   // run DFS and get Min-Cut
	void dfs(int s);
	void getMinCut(int s);
                   // return visited or not;
    int isVisited(int node_id);
};

}
