#pragma once

#include "load_inst.hpp"
#include "utility.hpp"
#include "mdd.hpp"


struct psa_nod {
	double val;
	size_t gap;
	double pes_val;
	size_t pes_gap;

	psa_nod(double _val, size_t _gap, double _pesv, size_t _pesg) { val = _val; gap = _gap; pes_val = _pesv; pes_gap = _pesg; }

};

void Build_mdd();

extern vector<vector<int> > value_ref;
extern vector<vector<double> > ext_gap_ref;
extern vector<vector<vector<psa_nod> >*> psa_mats;

extern vector<MDD*> MDDs;

extern double opt_psa;

