#pragma once

#include "node.hpp"
#include <ilcplex/ilocplex.h>

void build_model(IloModel& model);
void assoc_arcs();

extern IloEnv env;

extern bool SP_def, MP_def;

bool build_SP(IloRange& unb_ray);

class Var {
public:
	Arc* arc;
	IloNumVar var;
	double val;

	Var(Arc* _arc) { arc = _arc; var = IloNumVar(env, 0, 1, ILOFLOAT); val = 0; }

};

extern vector<Var> vars;


