#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <random>
#include <cassert>
#include <chrono>
#include "dem.h"

using namespace std;
using namespace std::chrono;

set<int> init(PROBLEM problem) {
	cout << "Starting initialization" << endl;
	set<int> init_sol;
	int best_i = -1;
	int best_score = 0;
	int fit = INFINITY;
	do {
		int best_i = -1; 
		int best_score = 0;
		int* edge_monitored_by = monitored_by(problem, init_sol);
		for (int i = 0; i < problem.nv; i++) {
			if (init_sol.find(i) != init_sol.end()) // it is already inside so skip it
				continue;
			int score = 0;
			for (int j = 0; j < problem.ne; j++)
				if (edge_monitored_by[j]==-1 && problem.e_on_any_path[i][j])
					score++;
			if (score > best_score) {
				best_score = score;
				best_i = i;
			}
		}
		if (best_i != -1)
			init_sol.insert(best_i);
		else
			break;
		int fit = fitness(problem, init_sol);
		cout << "Initial solution " << init_sol.size() <<" feasibility " << fit << " improvement "<< best_score << endl; // 
	} while (true);
	return init_sol;
}

// swapping k vertices and removing one
set<int> shake(set<int> sol, int k, int nv) {
	vector<int> shaked_sol;
	vector<int> shaked_sol_comp;
	for (int i = 0; i < nv; i++)
		if (sol.find(i) != sol.end())
			shaked_sol.push_back(i);
		else
			shaked_sol_comp.push_back(i);

	random_shuffle(shaked_sol.begin(), shaked_sol.end());
	random_shuffle(shaked_sol_comp.begin(), shaked_sol_comp.end());

	if (shaked_sol.size()>0 && shaked_sol_comp.size()>0) {
		// removes k vertices
		shaked_sol.erase(shaked_sol.end() - k, shaked_sol.end());
		// adds k vertices from solution complement
		for (int i = 0; i < k; i++)
			shaked_sol.push_back(shaked_sol_comp[i]);
	}
	else
		cout << "This shouldn't happen." << endl;

	// remove one random vertex -- to decrease candidate solution cardinality
	shaked_sol.erase(shaked_sol.begin() + rand() % shaked_sol.size());
	set<int> shaked_sol_set;
	for (int v : shaked_sol)
		shaked_sol_set.insert(v);
	return shaked_sol_set;
}

void vns(PROBLEM problem, int iter_max, int iter_max_wo_impr, int time_max, int k_min, int k_max_init, double p_move) {
	auto start = high_resolution_clock::now();

	calculate_edge_on_paths(problem);
	//print_problem(problem);

	set<int> best_sol = init(problem);
	int best_fit = fitness(problem, best_sol);
	// it must be feasible
	assert(best_fit == 0);

	int k = k_min;
	int last_impr_it = 0;
	int iter = 0;
	// vns main loop
	while(iter < iter_max && (iter-last_impr_it)<iter_max_wo_impr && duration_cast<seconds>(high_resolution_clock::now() - start).count()<time_max) {
		// adjusting k_max in each iteration -- // k_max can be larger than |sol| or |V/sol| 
		int k_max = min(min(k_max_init, (int)(best_sol.size())), int(problem.nv - best_sol.size()));
		set<int> new_sol = shake(best_sol, k, problem.nv);
		int fit = fitness(problem, new_sol);
		int new_fit = local_search_first_improvement(problem, new_sol);

		if (new_fit == 0) {
			// move to new solution
			best_sol = new_sol;
			best_fit = new_fit;
			k = k_min;
			last_impr_it = iter;
		}
		else {
			k++;
			if (k > k_max)
				k = k_min;
		}

		cout << "it=" << iter <<"\tt="<< duration_cast<seconds>(high_resolution_clock::now() - start).count()<< "s\tk=" << k << "\tk_max=" << k_max << "\tbest=" << best_sol.size() << endl;
		/*"\t{";
		for (int i : best_sol)
			cout << i << " ";
		cout << "}" << endl;*/
		iter++;
	}

	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - start);

	ofstream outf;
	outf.open("out_vns.txt", ios_base::app);
	outf << problem.name <<"\tnv=" << problem.nv << "\tne=" << problem.ne<< "\titer_max=" << iter_max << "\titer_wo_impr="<<iter_max_wo_impr<< "\tk_min=" << k_min << "\titer="<<iter<< "\tt=" << duration.count()/1000.0 << "\tk_max=" << k_max_init << "\tsize=" << best_sol.size() << "\tfeasible=" << (best_fit == 0) << "\tsol={";
	for (int i : best_sol)
		outf << i << " ";
	outf << "}" << endl;
	outf.close();
}

int main(int argc, char** argv)
{
	if (argc < 10) {
		cerr << "You must supply parameters <input_path> <starting_index> <iter_max> <iter_max_wo_impr> <time_max> <k_min> <k_max> <p_move> <seed>" << endl;
		exit(1);
	}

	char* fname = argv[1];
	int starting_index = atoi(argv[2]);
	int iter_max = atoi(argv[3]);
	int iter_max_wo_impr = atoi(argv[4]);
	int time_max = atoi(argv[5]);
	int k_min = atoi(argv[6]);
	int k_max = atoi(argv[7]);
	double p_move = atof(argv[8]);
	int seed = atoi(argv[9]);

	if (k_min < 1)
		throw invalid_argument("k_min should be at least 1.");
	if (k_min >= k_max)
		throw invalid_argument("k_min should be smaller than k_max");

	srand(seed);
	PROBLEM problem = load_problem(fname,starting_index);
	if (k_max > problem.nv) {
		cout << "Overriding k_max to " << problem.nv << endl;
		k_max = problem.nv;
	}

	vns(problem, iter_max, iter_max_wo_impr, time_max, k_min, k_max, p_move);

	return 0;
}
