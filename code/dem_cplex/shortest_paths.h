#pragma once
#include <vector>
using namespace std;

void add_edge(vector<vector<int>>& adj, int src, int dest);

vector<vector<int>> get_paths(vector<vector<int>> adj, int n, int start, int end);