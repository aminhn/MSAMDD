#pragma once

#include "load_inst.hpp"
#include "build_mdd.hpp"

bool filter(double opt_psa, double optimistic, double LP, int rew_pos, bool restart, bool reactive, int strt = 0);
void set_psa_del();

extern bool psa_feas, LP_feas;
