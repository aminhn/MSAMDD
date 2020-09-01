#include "build_model.hpp"
#include "build_mdd.hpp"

void build_expr(Node* nod, IloExpr& flow_cons);

vector<Var> vars;

bool SP_def = 0, MP_def = 0;

void build_model(IloModel& model) {

	assoc_arcs();

	IloExpr objexpr(env);

	for (vector<Var>::iterator it = vars.begin(); it != vars.end(); ++it) {
		if (it->arc->active)
			objexpr += it->arc->reward[0] * it->var;
	}

	model.add(IloMaximize(env, objexpr));

	objexpr.end();

	for (size_t i = 0; i < MDDs.size(); ++i) {
		bool first = 1;
		for (size_t j = 0; j < MDDs[i]->layers.size() - 1; ++j) {
			for (vector<vector<Node*>*>::iterator it2 = MDDs[i]->layers[j]->nodes.begin(); it2 != MDDs[i]->layers[j]->nodes.end(); ++it2) {
				for (vector<Node*>::iterator it = (*it2)->begin(); it != (*it2)->end(); ++it) {
					if ((*it)->active) {
						IloExpr flow_cons(env);

						build_expr((*it), flow_cons);

						if (first)
							model.add(flow_cons == 1);
						else
							model.add(flow_cons == 0);

						flow_cons.end();
					}
				}
			}
			first = 0;
		}

		/*IloExpr flow_cons(env);
		for (vector<Node*>::iterator it = MDDs[i]->layers.back()->nodes.begin(); it != MDDs[i]->layers.back()->nodes.end(); ++it) {
			if ((*it)->active) {
				vector<IloExpr> gap_cons;
				gap_cons.reserve((*it)->egap_children.size() + 1);

				build_expr((*it), flow_cons, gap_cons, vars, env);

				for (size_t k = 0; k < gap_cons.size(); ++k)
					model.add(gap_cons[k] <= 0);
			}
		}
		//cout << MDDs[i]->layers.size() << "\nOUT " << flow_cons << endl;
		//cin.get();

		model.add(flow_cons == -1);
		*/
	}

	//cout << model << endl;
	//cin.get();
}

void build_expr(Node* nod, IloExpr& flow_cons) {

	for (vector<Arc*>::iterator it2 = nod->children.begin(); it2 != nod->children.end(); ++it2) {
		if ((*it2)->active)
			flow_cons += vars[(*it2)->ID].var;
	}

	for (vector<Arc*>::iterator it2 = nod->parents.begin(); it2 != nod->parents.end(); ++it2) {
		if ((*it2)->active) 
			flow_cons -= vars[(*it2)->ID].var;
	}

	if (nod->gap_parent && nod->gap_parent->active) 
		flow_cons -= vars[nod->gap_parent->ID].var;
	
	if (nod->gap_child && nod->gap_child->active) 
		flow_cons += vars[nod->gap_child->ID].var;

}


void assoc_arcs() {

	cout << "Associating arcs\n";

	vars.clear();
	vars.shrink_to_fit();
	vars.reserve(num_act_arcs);

	size_t arc_count = 0;// nod_count = 0;
	for (size_t i = 0; i < MDDs.size(); ++i) {
		for (vector<Layer*>::iterator lay = MDDs[i]->layers.begin(); lay != MDDs[i]->layers.end(); ++lay) {
			for (vector<vector<Node*>*>::iterator it = (*lay)->nodes.begin(); it != (*lay)->nodes.end(); ++it) {
				for (vector<Node*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2) {
					if ((*it2)->active) {
						//(*it2)->ID = nod_count;
						//++nod_count;
						for (vector<Arc*>::iterator arc = (*it2)->children.begin(); arc != (*it2)->children.end(); ++arc) {
							if ((*arc)->active) {
								(*arc)->ID = arc_count;
								vars.emplace_back(*arc);
								++arc_count;
							}
						}
						if ((*it2)->gap_child && (*it2)->gap_child->active) {
							(*it2)->gap_child->ID = arc_count;
							vars.emplace_back((*it2)->gap_child);
							++arc_count;
						}
					}
				}
			}
		}
	}

	cout << "reality check: " << vars.size() << " " << num_act_arcs << endl;

}


bool build_SP(IloRange& unb_ray) {

	cout << "Building SP\n";

	IloEnv env_sp;
	IloModel model_sp(env_sp);

	vector<IloNumVarArray> col_vars;
	vector<IloRange> ord_cons, col_cons_p, col_cons_n;

	col_vars.reserve(strings.size());
	for (size_t i = 0; i < strings.size(); ++i) {
		IloNumVarArray var_row(env_sp, strings[i]->size(), 0, tot_chars + 1, ILOFLOAT);
		col_vars.emplace_back(var_row);
	}

	ord_cons.reserve(tot_chars - N);
	for (int i = 0; i < strings.size(); ++i) {
		for (int j = 0; j < strings[i]->size() - 1; ++j) {
			ord_cons.emplace_back(env_sp, col_vars[i][j] - col_vars[i][j + 1], -1);
			model_sp.add(ord_cons.back());
		}
	}

	col_cons_p.reserve(vars.size());
	col_cons_n.reserve(vars.size());
	for (vector<MDD*>::iterator it = MDDs.begin(); it != MDDs.end(); ++it) {
		for (int i = 1; i < (*it)->layers.size(); ++i) {
			for (int j = 0; j < (*it)->layers[i]->nodes[0]->size() - 1; ++j) {
				double sum = 0;
				for (vector<Arc*>::iterator it2 = (*it)->layers[i]->nodes[0]->at(j)->parents.begin(); it2 != (*it)->layers[i]->nodes[0]->at(j)->parents.end(); ++it2) {
					if ((*it2)->active)
						sum += vars[(*it2)->ID].val;
				}
				//cout << sum << " " << i << " " << j << " " << (*it)->layers.size() << endl;
				col_cons_p.emplace_back(env_sp, col_vars[(*it)->seq_fr][strings[(*it)->seq_fr]->size() - i] - col_vars[(*it)->seq_to][j], tot_chars * (1 - sum));
				col_cons_n.emplace_back(env_sp, -col_vars[(*it)->seq_fr][strings[(*it)->seq_fr]->size() - i] + col_vars[(*it)->seq_to][j], tot_chars * (1 - sum));
				model_sp.add(col_cons_p.back());
				model_sp.add(col_cons_n.back());
			}
		}
	}

	cout << "SP built\n";

	model_sp.add(IloMaximize(env_sp, 0));


	IloCplex SubP(model_sp);
	
	SubP.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	SubP.setParam(SubP.RootAlg, 1);
	SubP.setOut(env_sp.getNullStream());

	SubP.solve();
	double aa = 0, bb = 0;
	if (SubP.getStatus() != IloAlgorithm::Optimal) {
		cout << "Solution is infeasible\n";
		double LHS = 0;
		double RHS = 0;
		IloExpr expr(env);
		cout << "SIZE: " << ord_cons.size() << " " << col_cons_n.size() << endl;
		for (vector<IloRange>::iterator it = ord_cons.begin(); it != ord_cons.end(); ++it) {
			RHS += abs(SubP.getDual(*it));
		}

		size_t ii = 0;
		for (vector<MDD*>::iterator it = MDDs.begin(); it != MDDs.end(); ++it) {
			for (int i = 1; i < (*it)->layers.size(); ++i) {
				for (int j = 0; j < (*it)->layers[i]->nodes[0]->size() - 1; ++j) {
					double sum = 0;
					//cout << SubP.getDual(col_cons_p[ii]) << " " << SubP.getDual(col_cons_n[ii]) << endl;
					//cin.get();
					sum += tot_chars * -abs(SubP.getDual(col_cons_p[ii]));
					sum += tot_chars * -abs(SubP.getDual(col_cons_n[ii]));
					RHS += sum;
					for (vector<Arc*>::iterator it2 = (*it)->layers[i]->nodes[0]->at(j)->parents.begin(); it2 != (*it)->layers[i]->nodes[0]->at(j)->parents.end(); ++it2) {
						if ((*it2)->active) {
							expr += sum * vars[(*it2)->ID].var;
							aa += sum * vars[(*it2)->ID].val;
						}
					}
					//cout << sum << " " << i << " " << j << " " << (*it)->layers.size() << endl;
					++ii;
				}
			}
		}


		cout << endl << aa << " >= " <<   RHS + LHS << endl;

		unb_ray = IloRange(env, RHS + LHS, expr);

		//cout << unb_ray << endl;
		//cin.get();

		env_sp.end();

		return 0;

	}
	else {
		build_ord_mat(SubP, col_vars);
		//cin.get();
	}

	env_sp.end();

	return 1;

}


