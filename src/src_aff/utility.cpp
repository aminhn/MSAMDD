#include "utility.hpp"
#include "load_inst.hpp"
#include "build_mdd.hpp"
#include "build_model.hpp"
#include "filter.hpp"

float give_time(clock_t kk) {
	float ll = ((float)kk) / CLOCKS_PER_SEC;
	return ll;
}

 double get_value(vector<vector<int>*>& seqs) {

	double val = 0; 
	int gap_len = 0;
	for (int i = 0; i < N - 1; ++i) {
		for (int j = i + 1; j < N; ++j) {
			for (int k = 0; k < seqs[0]->size(); k++) {
				if (!(seqs[i]->at(k) == 45 && seqs[j]->at(k) == 45)) {
					if (seqs[i]->at(k) == 45 || seqs[j]->at(k) == 45)
						++gap_len;
					else if (gap_len == 0)
						val += subsmat[subsref[seqs[i]->at(k)]]->at(subsref[seqs[j]->at(k)]);
					else {
						val += subsmat[subsref[seqs[i]->at(k)]]->at(subsref[seqs[j]->at(k)]) - (opn_pen + ext_pen * gap_len);
						gap_len = 0;
					}
				}
			}
			if (gap_len > 0) {
				val -= opn_pen + ext_pen * gap_len;
				gap_len = 0;
			}
		}
	}

	return val;

}


void get_ext_ref(vector<vector<double> >& vec) {			//vector giving the extention cost of extending a gap of length "row number" by a gap of length "column number"
	vec.reserve(M + 2);
	vector<double> temp_vec;
	temp_vec.reserve(M + 2);
	vec.emplace_back(temp_vec);
	vec[0].emplace_back(0);
	for (int i = 1; i < M + 2; ++i) 
		vec[0].emplace_back(-(opn_pen + ext_pen * i));
	vector<double> temp_vec2;
	temp_vec2.reserve(M + 2);
	vec.emplace_back(temp_vec2);
	vec[1].emplace_back(0);
	for (int i = 1; i < M + 2; ++i)
		vec[1].emplace_back(-ext_pen * i);
}

void get_align_values(vector<int>& str2, vector<int>& val_ref, int pos) {
	val_ref.reserve(str2.size());
	for (int j = 0; j < str2.size(); j++) {
		val_ref.push_back(subsmat[pos]->at(subsref[str2[j]]));
	}
}

bool sol_exists(vector<vector<int>*>& mach_mat) {

	int mdd_count = 0;
	double sol_val = 0;
	for (int i = 0; i < mach_mat.size() - 1; ++i) {
		for (int j = i + 1; j < mach_mat.size(); ++j) {
			Node* cur_node = MDDs[mdd_count]->layers[0]->nodes[0]->at(0);
			bool gap = 0;
			int tnod = strings[j]->size();
			for (int s = mach_mat[i]->size() - 1; s >= 0; --s) {
				if (mach_mat[i]->at(s) != 45) {
					if (mach_mat[j]->at(s) == 45) {
						if ((!gap && (!cur_node->ogap_child || !cur_node->ogap_child->active)) || (gap && (!cur_node->egap_child || !cur_node->egap_child->active)))
							return 0;
						else {
							if (!gap) {
								gap = 1;
								sol_val += cur_node->ogap_child->reward[0];
								cur_node = cur_node->ogap_child->child;
							}
							else {
								sol_val += cur_node->egap_child->reward[0];
								cur_node = cur_node->egap_child->child;
							}
						}
					}
					else {
						if (cur_node->children.empty()) {
							return 0;
						}
						else {
							bool found = 0;
							for (vector<Arc*>::iterator it = cur_node->children.begin(); it != cur_node->children.end(); ++it) {
								if ((*it)->active && (*it)->child->state == tnod) {
									sol_val += (*it)->reward[0];
									//cout << (*it)->reward[0] << " " << (*it)->parent->state << " " << (*it)->child->state << "   ";
									cur_node = (*it)->child;
									--tnod;
									found = 1;
									gap = 0;
									break;
								}
							}
							if (!found)
								return 0;
						}
					}
				}
				else if (mach_mat[j]->at(s) != 45)
						--tnod;
			}
			++mdd_count;
			//cin.get();
		}
	}
	cout << "solution value in MDDs: " << sol_val << endl;

	return 1;
}

bool get_MIPstart(IloNumVarArray& mipvars, IloNumArray& mipvals) {

	int mdd_count = 0, strt = 0, end;
	for (int i = 0; i < heur_sol.size() - 1; ++i) {
		for (int j = i + 1; j < heur_sol.size(); ++j) {
			Node* cur_node = MDDs[mdd_count]->layers[0]->nodes[0]->at(0);
			bool gap = 0;
			int tnod = strings[j]->size();
			for (int s = heur_sol[i]->size() - 1; s >= 0; --s) {
				if (heur_sol[i]->at(s) != 45) {
					if (heur_sol[j]->at(s) == 45) {
						if ((!gap && (!cur_node->ogap_child || !cur_node->ogap_child->active)) || (gap && (!cur_node->egap_child || !cur_node->egap_child->active)))
							return 0;
						else {
							if (!gap) {
								gap = 1;
								end = cur_node->ogap_child->ID;
								for (int v = strt; v < end; ++v) {
									mipvars.add(vars[v].var);
									mipvals.add(0);
								}
								strt = end + 1;
								mipvars.add(vars[end].var);
								mipvals.add(1);
								cur_node = cur_node->ogap_child->child;
							}
							else {
								end = cur_node->egap_child->ID;
								for (int v = strt; v < end; ++v) {
									mipvars.add(vars[v].var);
									mipvals.add(0);
								}
								strt = end + 1;
								mipvars.add(vars[end].var);
								mipvals.add(1);
								cur_node = cur_node->egap_child->child;
							}
						}
					}
					else {
						if (cur_node->children.empty()) {
							return 0;
						}
						else {
							bool found = 0;
							for (vector<Arc*>::iterator it = cur_node->children.begin(); it != cur_node->children.end(); ++it) {
								if ((*it)->active && (*it)->child->state == tnod) {
									end = (*it)->ID;
									for (int v = strt; v < end; ++v) {
										mipvars.add(vars[v].var);
										mipvals.add(0);
									}
									strt = end + 1;
									mipvars.add(vars[end].var);
									mipvals.add(1);
									cur_node = (*it)->child;
									--tnod;
									found = 1;
									gap = 0;
									break;
								}
							}
							if (!found)
								return 0;
						}
					}
				}
				else if (heur_sol[j]->at(s) != 45)
					--tnod;
			}
			++mdd_count;
		}
	}

	return 1;
}


void get_viol(list<IloRange>& viol_cons) {
	
	vector<vector<int>> neig(comul_siz.back()), arc_ref(comul_siz.back());
	set<tuple<int, int, int>> cons_map;

	for (vector<Var>::iterator it = vars.begin(); it != vars.end(); ++it) {
		if (it->val > 0 && it->arc->parent->state != it->arc->child->state) {
			//cout << MDDs[it->arc->parent->mdd_count]->seq_fr << " " << MDDs[it->arc->parent->mdd_count]->seq_to << " " << strings[MDDs[it->arc->parent->mdd_count]->seq_fr]->size() - it->arc->child->lvl + 1 << " " << it->arc->child->state << endl;// << comul_siz[MDDs[it->arc->parent->mdd_count]->seq_to] + it->arc->child->state - 1<< endl;
			neig[comul_siz[MDDs[it->arc->parent->mdd_count]->seq_fr] + strings[MDDs[it->arc->parent->mdd_count]->seq_fr]->size() - it->arc->child->lvl].push_back(comul_siz[MDDs[it->arc->parent->mdd_count]->seq_to] + it->arc->child->state - 1);
			neig[comul_siz[MDDs[it->arc->parent->mdd_count]->seq_to] + it->arc->child->state - 1].push_back(comul_siz[MDDs[it->arc->parent->mdd_count]->seq_fr] + strings[MDDs[it->arc->parent->mdd_count]->seq_fr]->size() - it->arc->child->lvl);
		}
	}

	/*for (int i = 0; i < neig.size(); ++i) {
		cout << "N " << i << " : ";
		for (int j = 0; j < neig[i].size(); ++j) {
			cout << neig[i][j] << " ";
		}
		cout << endl;
	}*/
	
	for (int i = 0; i < neig.size(); ++i) {
		for (int j = 0; j < neig[i].size(); ++j) {
			for (int k = 0; k < neig[neig[i][j]].size(); k++) {
				if (neig[neig[i][j]][k] == i)
					continue;
				for (int p = 0; p < neig[i].size(); ++p) {
					if (neig[neig[i][j]][k] == neig[i][p])
						break;
					else if (nod_to_str[neig[neig[i][j]][k]] != nod_to_str[i] && (neig[neig[i][j]][k] < neig[i][p] || p == neig[i].size() - 1)) {
						if (i < neig[neig[i][j]][k])
							cons_map.emplace(neig[i][j], i, neig[neig[i][j]][k]);
						else
							cons_map.emplace(neig[i][j], neig[neig[i][j]][k], i);
						break;
					}
				}					
			}
		}
	}

	for (set<tuple<int, int, int>>::iterator it = cons_map.begin(); it != cons_map.end(); ++it) {
		int str1 = nod_to_str[get<0>(*it)], str2 = nod_to_str[get<1>(*it)], str3 = nod_to_str[get<2>(*it)];
		int state1 = get<0>(*it) - comul_siz[str1], state2 = get<1>(*it) - comul_siz[str2], state3 = get<2>(*it) - comul_siz[str3];
		
		Node* nod1;	Node* nod2;	Node* nod3;
		if (str1 < str2)
			nod1 = MDDs[comul_mdd[str1] + str2 - str1 - 1]->layers[strings[str1]->size() - state1]->nodes[0]->at(state2);
		else
			nod1 = MDDs[comul_mdd[str2] + str1 - str2 - 1]->layers[strings[str2]->size() - state2]->nodes[0]->at(state1);

		if (str1 < str3)
			nod2 = MDDs[comul_mdd[str1] + str3 - str1 - 1]->layers[strings[str1]->size() - state1]->nodes[0]->at(state3);
		else
			nod2 = MDDs[comul_mdd[str3] + str1 - str3 - 1]->layers[strings[str3]->size() - state3]->nodes[0]->at(state1);

		if (str2 < str3)
			nod3 = MDDs[comul_mdd[str2] + str3 - str2 - 1]->layers[strings[str2]->size() - state2]->nodes[0]->at(state3);
		else
			nod3 = MDDs[comul_mdd[str3] + str2 - str3 - 1]->layers[strings[str3]->size() - state3]->nodes[0]->at(state2);

		bool n1 = 0, n2 = 0;
		IloExpr expr(env);
		for (vector<Arc*>::iterator it2 = nod1->parents.begin(); it2 != nod1->parents.end(); ++it2) {
			if ((*it2)->active) {
				expr += vars[(*it2)->ID].var;
				n1 = 1;
			}
		}
		for (vector<Arc*>::iterator it2 = nod2->parents.begin(); it2 != nod2->parents.end(); ++it2) {
			if ((*it2)->active) {
				expr += vars[(*it2)->ID].var;
				n2 = 1;
			}
		}

		for (vector<Arc*>::iterator it2 = nod3->parents.begin(); it2 != nod3->parents.end(); ++it2) {
			if ((*it2)->active) {
				expr -= vars[(*it2)->ID].var;
			}
		}
		
		if (n1 && n2)
			viol_cons.emplace_back(env, expr, 1);
		expr.end();
	}
	
}

/*double run_MUSCLE(string& inst, bool refine) {

	string str;
	if (refine)
		str = "cd " + path_file + " & muscle.exe -in " + inst + ".fa -out " + inst + "_s.txt -matrix blosum.ncbi -gapopen " + to_string(-opn_pen) + " -gapextend " + to_string(-ext_pen) + " -center 0.0 -sp  -refine -quiet";
	else
		str = "cd " + path_file + " & muscle.exe -in " + inst + ".fa -out " + inst + "_s.txt -matrix blosum.ncbi -gapopen " + to_string(-opn_pen) + " -gapextend " + to_string(-ext_pen) + " -center 0.0 -sp -quiet";
	
	const char *command = str.c_str();
	bool exc = system(command);

	if (exc)
		cout << "MUSCLE did not execute! Using previous solution if available\n";

	string out_f = inst + "_s.txt";

	return Load_sol(out_f);

}*/

double run_MUSCLE(string& inst, bool refine) {

	string str;
	if (refine)
		str = "chmod +x ./muscle && ./muscle -in " + inst + ".fa -out " + inst + "_s.txt -matrix " + subs_file + " -gapopen " + to_string(-opn_pen-ext_pen) + " -gapextend " + to_string(-ext_pen) + " -center 0.0 -refine -quiet";
	else
		str = "chmod +x ./muscle && ./muscle -in " + inst + ".fa -out " + inst + "_s.txt -matrix " + subs_file + " -gapopen " + to_string(-opn_pen-ext_pen) + " -gapextend " + to_string(-ext_pen) + " -center 0.0 -quiet";

	const char *command = str.c_str();
	bool exc = system(command);

	if (exc)
		cout << "MUSCLE did not execute! Using previous solution if available\n";

	string out_f = inst + "_s.txt";

	return Load_sol(out_f);

}


void build_ord_mat(IloCplex& SubP, vector<IloNumVarArray>& col_vars) {
	
	vector<vector<int>> ord_mat;
	ord_mat.reserve(N);

	size_t j = 0;
	for (vector<IloNumVarArray>::iterator it = col_vars.begin(); it != col_vars.end(); ++it) {
		ord_mat.push_back(vector<int>(tot_chars, 45));
		for (int i = 0; i < it->getSize(); ++i) {
			ord_mat.back()[SubP.getValue((*it)[i])] = strings[j]->at(i);
		}
		++j;
	}

	vector<bool> indic(tot_chars, 0);
	for (int i = 0; i < ord_mat[0].size(); ++i) {
		for (int j = 0; j < ord_mat.size(); ++j) {
			if (ord_mat[j][i] != 45) {
				indic[i] = 1;
				break;
			}
		}
	}

	vector<vector<int>*> output;
	for (int i = 0; i < ord_mat.size(); ++i) {
		output.emplace_back(new vector<int>);
		for (int j = 0; j < ord_mat[i].size(); ++j) {
			if (indic[j])
				output[i]->push_back(ord_mat[i][j]);
		}
	}

	double sol_val = get_value(output);

	cout << "Solution value is: " << sol_val << endl;

	if (sol_val > incumb) {
		disp_sol(output);
		check_feas(output);
		cout << "redefining heur mat\n";
		for (int i = 0; i < heur_sol.size(); ++i) {
			delete heur_sol[i];
			heur_sol[i] = output[i];
		}
		incumb = sol_val;
		cout << "Incumbent updated to new found solution value!\n";
		for (vector<MDD*>::iterator it = MDDs.begin(); it != MDDs.end(); ++it) {
			(*it)->comp_lng_to(0, 0);
			(*it)->comp_lng_fr(0, 0);
		}
		filter(opt_psa, incumb, 0, 0, 1, 0);
		
		psa_feas = sol_exists(output);
		LP_feas = psa_feas;

		cout << "solution is in " << psa_feas << endl;

	}
	else {
		for (int i = 0; i < output.size(); ++i)
			delete output[i];
	}

}

void check_feas(vector<vector<int>*> mach_mat) {

	for (size_t i = 0; i < N; ++i) {
		size_t pos = 0;
		for (size_t j = 0; j < mach_mat[0]->size(); ++j) {
			if (mach_mat[i]->at(j) == 45)
				continue;
			else if (mach_mat[i]->at(j) != strings[i]->at(pos)) {
				cout << "NOT OK!!!!!!\n";
				break;
			}
			++pos;
		}
		if (pos == strings[i]->size())
			cout << "OK\n";
	}
}

void disp_sol(vector<vector<int>*> sol) {

	for (int i = 0; i < sol.size(); ++i) {
		for (int j = 0; j < sol[i]->size(); ++j) {
			cout << (char)sol[i]->at(j);
		}
		cout << endl;
	}
}
