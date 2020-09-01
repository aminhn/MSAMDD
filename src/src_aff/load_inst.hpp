#pragma once

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <time.h>
#include <map>
#include <limits>
#include <math.h>
#include <list>

using namespace std;

bool Load_instance(string&);
bool Load_subsmat();
double Load_sol(string&);
void write_sol(list<vector<int> >&, string&);

extern vector<vector<int>*> strings, subsmat, heur_sol;
extern vector<int> alphabet, nod_to_str, comul_siz, comul_mdd, num_gap_arcs, time_limit;
extern vector<double> elapsed_time, filt;
extern map<int, int> subsref;
extern vector<string> seq_names;


extern size_t N, M, tot_chars;
extern int num_arcs, num_act_arcs, phaze;
extern double opn_pen, ext_pen, incumb;
extern string item_file, subs_file;

extern clock_t start_time;
