#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>

#include "shortest_paths.h"
#include "dem.h"
#include <algorithm>
#include <cassert>

using namespace std;


void print_problem(PROBLEM problem) {
	cout << problem.name;
	cout << problem.nv << "\t" << problem.ne << endl;
	// print adjacency matrix
	for (int i = 0; i < problem.nv; i++) {
		for (int j = 0; j < problem.nv; j++)
			cout << problem.adj_mat[i][j] << "\t";
		cout << endl;
	}
	// for each pair (v, e) print if there is u such that e is on all shortest paths between v and u
	for (int i = 0; i < problem.nv; i++) {
		cout << i << "\t";
		for (int k = 0; k < problem.ne; k++)
			if (problem.e_on_any_path[i][k])
				cout << "e" << k << "={" << problem.edges[k].first << "\t" << problem.edges[k].second << "}" << "\t";
		cout << endl;
	}
}

void calculate_edge_on_paths_slow(PROBLEM& problem) {
	cout << "Calculating shortest paths." << endl;
	problem.e_on_any_path = new bool* [problem.nv];
	for (int i = 0; i < problem.nv; i++) {
		problem.e_on_any_path[i] = new bool[problem.ne];
		for (int p = 0; p < problem.ne; p++)
			problem.e_on_any_path[i][p] = false;
	}
	// vector of vectors is used to store the graph in the form of an adjacency list
	vector<vector<int>> adj_list(problem.nv, vector<int>());
	// creating adjacency list based on the adjacency matrix
	for (int i = 0; i < problem.nv; i++)
		for (int j = 0; j < i; j++)
			if (problem.adj_mat[i][j])
				add_edge(adj_list, i, j);
	// now filling the 2D structure e_on_any_path for each v and each e
	vector<vector<int>> paths;
	vector<vector<int>> edge_counts(problem.nv, vector<int>(problem.nv));
	for (int i = 0; i < problem.nv; i++) {
		cout << "Paths from vertex " << i << endl;
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < problem.nv; k++)
				for (int r = 0; r < problem.nv; r++)
					edge_counts[k][r] = 0;
			paths = get_paths(adj_list, problem.nv, i, j);
			for (auto v : paths) {
				reverse(v.begin(), v.end());
				int prev = -1;
				for (int u : v) {
					if (prev != -1) {
						edge_counts[prev][u] += 1;
						edge_counts[u][prev] += 1;
					}
					prev = u;
				}
			}
			// now checking if some edge was on all shortest paths between i and some j or between j and some i
			for (int k = 0; k < problem.nv; k++)
				for (int r = 0; r < problem.nv; r++)
					if (edge_counts[k][r] == paths.size()) {
						// now find which edge is this in the edge order 
						for (int p = 0; p < problem.ne; p++)
							if ((problem.edges[p].first == k && problem.edges[p].second == r) || (problem.edges[p].first == r && problem.edges[p].second == k)) {
								problem.e_on_any_path[i][p] = true;
								problem.e_on_any_path[j][p] = true;
							}
					}
		}
	}
	cout << "Shortest path calculations finished." << endl;
}

// calculation explained in email 5.6.2023 7.37PM and Lemma 5 in Foucaud 2022. 
// Lemma 5. Let x be a vertex of a connected graph G. Then, an edge uv belongs to EM(x) if and only if u \in Li(x) and v is the only neighbour of u in Li−1(x), for some integer i.
// where Li(x) is the set of all vertices at distance i from vertex x. 
void calculate_edge_on_paths(PROBLEM& problem) {
	cout << "Calculating shortest distances with Floyd-Warshall." << endl;
	int** dist = new int* [problem.nv];
	for (int i = 0; i < problem.nv; i++)
		dist[i] = new int[problem.nv];
	for (int i = 0; i < problem.nv; i++)
		for (int j = 0; j < problem.nv; j++)
			if (i == j)
				dist[i][j] = 0;
			else if (problem.adj_mat[i][j])
				dist[i][j] = 1;
			else
				dist[i][j] = problem.nv+1;
	for (int k = 0; k < problem.nv; k++)
		for (int i = 0; i < problem.nv; i++)
			for (int j = 0; j < problem.nv; j++)
				if (dist[i][j] > dist[i][k] + dist[k][j])
					dist[i][j] = dist[i][k] + dist[k][j];

	//for (int i = 0; i < problem.nv; i++) {
	//	for (int j = 0; j < problem.nv; j++)
	//		cout << dist[i][j] << "\t";
	//	cout << endl;
	//}

	problem.e_on_any_path = new bool* [problem.nv];
	for (int i = 0; i < problem.nv; i++) {
		problem.e_on_any_path[i] = new bool[problem.ne];
		for (int p = 0; p < problem.ne; p++)
			problem.e_on_any_path[i][p] = false;
	}

	for (int e = 0; e < problem.ne; e++) {
		cout << "Checking if edge " << e << endl;
		for (int x = 0; x < problem.nv; x++) {
			// does edge e belongs to EM(x)
			int u = problem.edges[e].first;
			int v = problem.edges[e].second;
			if (dist[x][u] == dist[x][v])
				continue;	// node x cannot monitor edge e because both of e ending points are at the same distance -- if the optimal path goes through this edge this is never the case
			if (dist[x][u] == dist[x][v] + 1) {
				// checking if v is the only neighbor of u at distance dist[x][u]-1 from x
				bool maybe_x = true;
				for (int w = 0; w < problem.nv; w++) {
					if (w == v)
						continue;
					if (dist[x][u] == dist[x][w] + 1) {
						maybe_x = false; // x surely does not monitor e=vu because there is also an edge wu that is a part of some other shortest path
						break;
					}
				}
				if (!maybe_x)
					continue;
			}
			if (dist[x][v] == dist[x][u] + 1) {
				// checking if u is the only neighbor of v at distance dist[x][v]-1 from x
				bool maybe_x = true;
				for (int w = 0; w < problem.nv; w++) {
					if (w==u)
						continue;
					if (dist[x][v] == dist[x][w] + 1) {
						maybe_x = false; // x surely does not monitor e=uv because there is also an edge wv that is a part of some other shortest path
						break;
					}
				}
				if (!maybe_x)
					continue;
			}
			// x monitors edge e
			problem.e_on_any_path[x][e] = true;
		}
	}
	cout << "Shortest path calculations finished." << endl;
}

PROBLEM load_problem(char* fname, int starting_index) {
	PROBLEM problem;
	problem.name = fname;
	ifstream fs;
	fs.open(fname);
	fs >> problem.nv;
	fs >> problem.ne;
	problem.adj_mat = new bool* [problem.nv];
	for (int i = 0; i < problem.nv; i++) {
		problem.adj_mat[i] = new bool[problem.nv];
		for (int j = 0; j < problem.nv; j++)
			problem.adj_mat[i][j] = false;
	}

	int p, q;
	// each line contains edge endpoints
	for (int i = 0; i < problem.ne; i++) {
		fs >> p;
		fs >> q;
		p -= starting_index;
		q -= starting_index;
		problem.adj_mat[p][q] = true;
		problem.adj_mat[q][p] = true;
		if (p < q)
			problem.edges.push_back(pair<int, int>(p, q));
		else if (p > q)
			problem.edges.push_back(pair<int, int>(q, p));
		else
			throw invalid_argument("Self-loops are not allowed.");
	}

	fs.close();

	return problem;
}

int* monitored_by(PROBLEM problem, set<int> sol) {
	int* edge_monitored_by = new int[problem.ne];
	for (int i = 0; i < problem.ne; i++)
		edge_monitored_by[i] = -1;
	// now going through the monitoring set (sol) and checking what edges are monitored
	for (int k = 0; k < problem.ne; k++) {
		for (int i : sol) {
			if (problem.e_on_any_path[i][k]) {
				edge_monitored_by[k] = i;
				break;
			}
		}
	}
	return edge_monitored_by;
}

int fitness(PROBLEM problem, set<int> sol) {
	int* edge_monitored_by = monitored_by(problem, sol);
	// countin the number of monitored edges
	int monitored = 0;
	for (int i = 0; i < problem.ne; i++)
		if (edge_monitored_by[i] != -1)
			monitored += 1;
	int infeasibility = problem.ne - monitored; // when this is 0, the solution is feasible
	return infeasibility;
}

int local_search_first_improvement(PROBLEM problem, set<int>& sol) {
	int best_fit = fitness(problem, sol);
	bool impr = true;

	while (impr) {
		impr = false;
		set<int> new_sol(sol);
		vector<int> sol_vect;
		vector<int> sol_vect_comp;
		for (int i = 0; i < problem.nv; i++)
			if (sol.find(i) == sol.end())
				sol_vect_comp.push_back(i);
			else
				sol_vect.push_back(i);

		random_shuffle(sol_vect.begin(), sol_vect.end());
		random_shuffle(sol_vect_comp.begin(), sol_vect_comp.end());

		int* edge_monitored_by = monitored_by(problem, sol);
		vector<vector<int>> edges_monitored_by_vertex(problem.nv, vector<int>());
		vector<int> edges_not_monitored;
		for (int k = 0; k < problem.ne; k++)
			if (edge_monitored_by[k] == -1)
				edges_not_monitored.push_back(k);
			else
				edges_monitored_by_vertex[edge_monitored_by[k]].push_back(k);

		if (edges_not_monitored.size() == 0) // there is no chance for improvement when everything is already covered, diff cannot become negative
			break;

		// 1-swap first improvement 
		for (int i : sol_vect) {
			for (int j : sol_vect_comp) {
				// calling method and passing too many arguments data slows down a program significantly
				//int fit_diff =  fitness_difference(problem, sol, edges_monitored_by_vertex[i], edges_not_monitored, i, j);
				// checking only edges that are either not monitored or that are monitored with a node that is currently being removed
				int diff = 0;
				for (auto e : edges_not_monitored)
					if (problem.e_on_any_path[j][e])
						diff--;

				for (auto e : edges_monitored_by_vertex[i]) {
					// it could have been monitored with some other vertex and not only newly added_idx so we check all other included vertices plus added_idx newly added minus removed_idx
					diff++; // first we assume it is not monitored anymore
					if (problem.e_on_any_path[j][e])
						diff--;	// this neutralizes the above diff++
					else {
						for (int s : sol) {
							if (s != i && problem.e_on_any_path[s][e]) {
								diff--; // it is monitored (neutralized) with somebody else in the solution set
								break;
							}
						}
					}
				}

				if (diff < 0) {
					sol.erase(i);
					sol.insert(j);
					best_fit+=diff;
					int test_fit = fitness(problem, sol);
					assert(test_fit == best_fit);
					impr = true;
					goto next_iter;
				}
			}
		}
	next_iter:;
	}
	return best_fit;
}

int local_search_best_improvement(PROBLEM problem, set<int>& sol) {
	int best_fit = fitness(problem, sol);
	int best_i, best_j;
	bool impr = true;

	while (impr) {
		impr = false;
		best_i = -1;
		best_j = -1;
		set<int> new_sol(sol);
		vector<int> sol_vect;
		vector<int> sol_vect_comp;
		for (int i = 0; i < problem.nv; i++)
			if (sol.find(i) == sol.end())
				sol_vect_comp.push_back(i);
			else
				sol_vect.push_back(i);

		// 1-swap best improvement 
		for (int i : sol_vect) {
			new_sol.erase(i);
			for (int j : sol_vect_comp) {
				new_sol.insert(j);
				int new_fit = fitness(problem, new_sol);
				if (new_fit < best_fit) {
					best_fit = new_fit;
					best_i = i;
					best_j = j;
					impr = true;
				}
				new_sol.erase(j);
			}
			new_sol.insert(i);
		}
		if (impr) {
			sol.erase(best_i);
			sol.insert(best_j);
			int test_fit = fitness(problem, sol);
			assert(test_fit == best_fit);
		}
	}
	return best_fit;
}