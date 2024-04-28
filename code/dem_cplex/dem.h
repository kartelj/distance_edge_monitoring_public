#pragma once

using namespace std;

typedef struct {
	string name;
	int nv; // number of vertices
	int ne; // number of edges
	bool** adj_mat; // adjacency matrix
	vector<pair<int, int>> edges;
	bool** e_on_any_path; // e_on_path[u][e] = true if e is on all shortest paths between nodes u and at least one node v from V(G)
} PROBLEM;

void print_problem(PROBLEM problem);

void calculate_edge_on_paths(PROBLEM& problem);

int* monitored_by(PROBLEM problem, set<int> sol);

PROBLEM load_problem(char* fname, int starting_index);

int fitness(PROBLEM problem, set<int> sol);

int local_search_first_improvement(PROBLEM problem, set<int>& sol);

int local_search_best_improvement(PROBLEM problem, set<int>& sol);