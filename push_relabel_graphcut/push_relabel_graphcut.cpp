#include "push_relabel_graphcut.hpp"

using namespace std;

namespace PushRelabel
{

Graph::Graph(int V)
{
	this->V = V;

	// all vertices are initialized with 0 height 
	// and 0 excess flow 
	for (int i = 0; i < V; i++)
		ver.push_back(Vertex(0, 0));

	std::vector<int> tmp;
	edge_map.assign(V, tmp);

	visited.assign(V, false);
}

Graph::Graph()
{
	this->V = 0;
}

void Graph::addEdge(int u, int v, int capacity)
{
	// flow is initialized with 0 for all edge 
	edge.push_back(Edge(0, capacity, u, v));
	edge_map[u].push_back(edge.size()-1);
}

int Graph::addNode()
{
    this->V += 1;
	ver.push_back(Vertex(0, 0));
    visited.push_back(false);

	std::vector<int> tmp;
	edge_map.push_back(tmp);

    return (this->V - 1);
}

void Graph::preflow(int s)
{
	// Making h of source Vertex equal to no. of vertices 
	// Height of other vertices is 0. 
	ver[s].h = ver.size();

	// Fast pointer access
	PushRelabel::Edge* edge_ptr = &edge[0];
	PushRelabel::Vertex* ver_ptr = &ver[0];
	int* edge_map_ptr = &edge_map[s][0];

	// If current edge goes from source 
	for(int i = 0; i < edge_map[s].size(); i++)
	{
		int idx = edge_map_ptr[i];
		// Flow is equal to capacity 
		edge_ptr[idx].flow = edge_ptr[idx].capacity;

		// Initialize excess flow for adjacent v 
		ver_ptr[edge_ptr[idx].v].e_flow += edge_ptr[idx].flow;

		// Add an edge from v to s in residual graph with 
		// capacity equal to 0 
		edge.push_back(Edge(-edge_ptr[idx].flow, 0, edge_ptr[idx].v, s));
		edge_map[edge_ptr[idx].v].push_back(edge.size()-1);
	}
}

// Update reverse flow for flow added on ith Edge 
void Graph::updateReverseEdgeFlow(int i, int flow)
{
	int u = edge[i].v, v = edge[i].u;

	// Fast pointer access
	PushRelabel::Edge* edge_ptr = &edge[0];
	int* edge_map_ptr = &edge_map[u][0];
	for (int j = 0; j < edge_map[u].size(); j++)
	{
		int idx = edge_map_ptr[j];
		if(edge_ptr[idx].v == v)
		{
			edge_ptr[idx].flow -= flow;
			return;
		}
	}

	// adding reverse Edge in residual graph 
	Edge e = Edge(0, flow, u, v);
	edge.push_back(e);
	edge_map[u].push_back(edge.size()-1);
}

// To push flow from overflowing vertex u 
bool Graph::push(int u)
{
	// Fast pointer access
	PushRelabel::Edge* edge_ptr = &edge[0];
	PushRelabel::Vertex* ver_ptr = &ver[0];
	int* edge_map_ptr = &edge_map[u][0];
	
	// Traverse hashmap of edge to find an adjacent (of u) 
	// to which flow can be pushed 
	// overflowing vertex 
	for (int i = 0; i < edge_map[u].size(); i++)
	{
		int idx = edge_map_ptr[i];

		// if flow is equal to capacity then no push 
		// is possible 
		if (edge_ptr[idx].flow == edge_ptr[idx].capacity)
			continue;

		// Push is only possible if height of adjacent 
		// is smaller than height of overflowing vertex 
		if (ver_ptr[u].h > ver_ptr[edge_ptr[idx].v].h)
		{
			// Flow to be pushed is equal to minimum of 
			// remaining flow on edge and excess flow. 
			int flow = std::min(edge_ptr[idx].capacity - edge_ptr[idx].flow, ver_ptr[u].e_flow);

			// Reduce excess flow for overflowing vertex 
			ver_ptr[u].e_flow -= flow;

			// Increase excess flow for adjacent 
			ver_ptr[edge_ptr[idx].v].e_flow += flow;

			// Add residual flow (With capacity 0 and negative 
			// flow) 
			edge_ptr[idx].flow += flow;

			updateReverseEdgeFlow(idx, flow);

			return true;
		}
	}

	return false;
}

// function to relabel vertex u 
void Graph::relabel(int u)
{
	// Initialize minimum height of an adjacent 
	int mh = INT_MAX;

	// Fast pointer access
	PushRelabel::Edge* edge_ptr = &edge[0];
	PushRelabel::Vertex* ver_ptr = &ver[0];
	int* edge_map_ptr = &edge_map[u][0];

	// Find the adjacent with minimum height 
	for (int i = 0; i < edge_map[u].size(); i++)
	{
		int idx = edge_map_ptr[i];
		// if flow is equal to capacity then no 
		// relabeling 
		if (edge_ptr[idx].flow == edge_ptr[idx].capacity)
			continue;

		// Update minimum height 
		if (ver_ptr[edge_ptr[idx].v].h < mh)
		{
			mh = ver_ptr[edge_ptr[idx].v].h;

			// updating height of u 
			ver_ptr[u].h = mh + 1;
		}
	}
}

// main function for printing maximum flow of graph 
int Graph::getMaxFlow(int s, int t)
{
	this->s = s;
	this->t = t;
	
	preflow(s);

	bool overflow_exist = true;
	while(overflow_exist)
	{
		overflow_exist = false;
		for(int i = 0; i < V; i++)
		{
			if((ver[i].e_flow > 0) && (i != s) && (i != t))
			{
				if (!push(i))
					relabel(i);
				overflow_exist = true;
			}
		}
	}

	// ver.back() returns last Vertex, whose 
	// e_flow will be final maximum flow 
	return ver[t].e_flow;
}

// A DFS based function to find all reachable vertices from s.  The function 
// marks visited[i] as true if i is reachable from s.  The initial values in 
// visited[] must be false. We can also use BFS to find reachable vertices 
void Graph::dfs(int s)
{
	visited[s] = true;

	// Fast pointer access
	PushRelabel::Edge* edge_ptr = &edge[0];
	int* edge_map_ptr = &edge_map[s][0];

	for (int i = 0; i < edge_map[s].size(); i++)
	{
		int idx = edge_map_ptr[i];
		if((edge_ptr[idx].capacity > 0) 
			&& !(edge_ptr[idx].capacity == edge_ptr[idx].flow) 
			&& !visited[edge_ptr[idx].v])
		{
			dfs(edge_ptr[idx].v);
		}
	}
}

void Graph::getMinCut(int s)
{
    // run dfs
	dfs(s);
}

int Graph::isVisited(int node_id)
{
    if(visited[node_id] == true)
        return 1;
    else
        return 0;
}

}