#pragma once
#include <string>
#include <vector>
#include "build_mdd.hpp"
#include <time.h>
#include <set>
#include <tuple>
#include <algorithm>
#include "build_model.hpp"


using namespace std;

float give_time(clock_t kk);
double get_value(vector<vector<int>*>&);
bool sol_exists(vector<vector<int>*>& mach_mat);
void get_align_values(vector<int>& str2, vector<int>& val_ref, int chr);
void get_ext_ref(vector<vector<double>>& vec);
void get_viol(list<IloRange>& viol_cons);
void build_ord_mat(IloCplex& SubP, vector<IloNumVarArray>& col_vars);
void check_feas(vector<vector<int>*> mach_mat);
void disp_sol(vector<vector<int>*> sol);
bool get_MIPstart(IloNumVarArray& mipvars, IloNumArray& mipvals);
double run_MUSCLE(string& inst, bool refine);


