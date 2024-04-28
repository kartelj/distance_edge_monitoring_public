#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <random>
#include <cassert>
#include <chrono>
#include "dem.h"
#include <ilcplex/ilocplexi.h>

using namespace std;
using namespace std::chrono;

void cplex_solve(PROBLEM problem, int time_max, int thread_cnt) {
	auto start = high_resolution_clock::now();

	calculate_edge_on_paths(problem);
	//print_problem(problem);

	IloEnv env;
	// cost matrix -- |V| x |E| binary matrix C[v, e] = true if exists u such that e belongs to all shortest paths between v and u
	IloArray<IloIntArray> C(env);
	for (int i = 0; i < problem.nv; i++) {
		IloIntArray arr = IloIntArray(env, problem.ne);
		for (int j = 0; j < problem.ne; j++)
			if (problem.e_on_any_path[i][j])
				arr[j]=1;
			else
				arr[j]=0;
		C.add(arr);
	}
	// variable if vertex is included inside the DEM
	IloNumVarArray dem(env, problem.nv, 0, 1, ILOINT);

	// Create the model.
	IloModel model(env);
	// Constraint: Each edge must be covered with at least one v from DEM such that there is an u for which e lies on all shortest paths between v and u
	for (IloInt i = 0; i < problem.ne; ++i){
		IloExpr vC(env);
		for (IloInt j = 0; j < problem.nv; j++)
			vC += dem[j]*C[j][i];
		model.add(vC >= 1);
		vC.end();
	}

	// Objective: Minimize the number of selected vertices
	IloExpr obj = IloSum(dem);
	model.add(IloMinimize(env, obj));
	obj.end();

	// Create a solver instance and extract the model to it.
	IloCplex cplex(env);
	cplex.setParam(IloCplex::Param::TimeLimit, time_max);
	cplex.setParam(IloCplex::Param::Threads, thread_cnt);
	cplex.extract(model);
	cplex.solve();

	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - start);

	ofstream outf;
	outf.open("out_cplex.txt", ios_base::app);
	outf << problem.name << "\tnv="<<problem.nv << "\tne="<<problem.ne<< "\ttime_max=" << time_max << "\tthread_cnt="<<thread_cnt<< "\tt=" << duration.count() / 1000.0 << "\tsize=" << cplex.getObjValue() << "\tstatus=" << cplex.getStatus() << "\tsol={";
	for (IloInt i=0; i<problem.nv; i++)
		if(cplex.getValue(dem[i])==1)
			outf << i << " ";
	outf << "}" << endl;
	outf.close();
}

int main(int argc, char** argv)
{
	if (argc < 5) {
		cerr << "You must supply parameters <input_path> <starting_index> <time_max> <thread_count>" << endl;
		exit(1);
	}

	char* fname = argv[1];
	int starting_index = atoi(argv[2]);
	int time_max = atoi(argv[3]);
	int thread_cnt = atoi(argv[4]);

	PROBLEM problem = load_problem(fname, starting_index);
	cplex_solve(problem, time_max, thread_cnt);

	return 0;
}
