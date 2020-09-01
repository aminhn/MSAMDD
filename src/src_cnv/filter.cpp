# include "filter.hpp"

int num_act_nodes;

bool filter(double opt_psa, double optimistic, double LP, int rew_pos, bool restart, bool reactive, int strt) {

	if (rew_pos == 0)
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&Performing PSA filtering\n";
	else
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&Performing LP filtering\n";

	bool init = 0;
	if (MDDs.size() - strt == 1 || optimistic == incumb)
		init = 1;

	int init_num = num_act_arcs;

	cout << "initial number of arcs: " << num_arcs << " number active: " << num_act_arcs << endl;

	vector<bool> del(vars.size(), 0);
	IloNumVarArray del_vars(env);

	double bound = LP;
	int num_inact_arc, num_inact_nod, num_inact_layer;
	num_act_nodes = 0;
	for (int s = strt; s < MDDs.size(); ++s) {
		//cout << "MDD " << s << endl;
		for (int i = MDDs[s]->layers.size() - 1; i >= 0; --i) {
			//cout << "layer " << i << endl;
			if (rew_pos == 0)
				bound = opt_psa - MDDs[s]->lng_path->max_to[0];

			num_inact_layer = 0;
			for (vector<vector<Node*>*>::iterator it = MDDs[s]->layers[i]->nodes.begin(); it != MDDs[s]->layers[i]->nodes.end(); ++it){
				num_inact_nod = 0;
				for (int n = (*it)->size() - 1; n >= 0; --n) {
					num_inact_arc = 0;
					Node* cur_nod = (*it)->at(n);

					for (int j = cur_nod->parents.size() - 1; j >= 0; --j) {
						if (restart || cur_nod->parents[j]->active) {
							if (init && cur_nod->parents[j]->reward[rew_pos] + cur_nod->parents[j]->parent->max_to[rew_pos] + cur_nod->max_fr[rew_pos] + 0.00001 < incumb - bound) {
								if (cur_nod->parents[j]->active) {
									--num_act_arcs;
									if (!vars.empty()) {
										del_vars.add(vars[cur_nod->parents[j]->ID].var);
										del[cur_nod->parents[j]->ID] = 1;
									}
								}
								cur_nod->parents.erase(cur_nod->parents.begin() + j);
								--num_arcs;
							}
							else if (cur_nod->parents[j]->reward[rew_pos] + cur_nod->parents[j]->parent->max_to[rew_pos] + cur_nod->max_fr[rew_pos] + 0.00001 < optimistic - bound) {
								++num_inact_arc;
								if (cur_nod->parents[j]->active) {
									cur_nod->parents[j]->active = 0;
									--num_act_arcs;
									if (rew_pos == 0)
										cur_nod->parents[j]->del_psa = 1;
									else {
										cur_nod->parents[j]->del_psa = 0;
									}
									if (!vars.empty()) {
										del_vars.add(vars[cur_nod->parents[j]->ID].var);
										del[cur_nod->parents[j]->ID] = 1;
									}
								}
							}
							else if (restart && !cur_nod->parents[j]->active && reactive && (rew_pos == 0 || !cur_nod->parents[j]->del_psa)) {
								cur_nod->parents[j]->active = 1;
								cur_nod->parents[j]->del_psa = 0;
								++num_act_arcs;
							}
						}
						else
							++num_inact_arc;
					}

					if (cur_nod->gap_parent && (restart || cur_nod->gap_parent->active)) {
						if (init && cur_nod->gap_parent->reward[rew_pos] + cur_nod->gap_parent->parent->max_to[rew_pos] + cur_nod->max_fr[rew_pos] + 0.00001 < incumb - bound) {
							if (cur_nod->gap_parent->active) {
								--num_act_arcs;
								if (!vars.empty()) {
									del_vars.add(vars[cur_nod->gap_parent->ID].var);
									del[cur_nod->gap_parent->ID] = 1;
								}
							}
							cur_nod->gap_parent->parent->gap_child = NULL;
							delete cur_nod->gap_parent;
							cur_nod->gap_parent = NULL;
							--num_arcs;
						}
						else if (cur_nod->gap_parent->reward[rew_pos] + cur_nod->gap_parent->parent->max_to[rew_pos] + cur_nod->max_fr[rew_pos] + 0.00001 < optimistic - bound) {
							if (cur_nod->gap_parent->active) {
								cur_nod->gap_parent->active = 0;
								--num_act_arcs;
								if (rew_pos == 0)
									cur_nod->gap_parent->del_psa = 1;
								else {
									cur_nod->gap_parent->del_psa = 0;
								}
								if (!vars.empty()) {
									del_vars.add(vars[cur_nod->gap_parent->ID].var);
									del[cur_nod->gap_parent->ID] = 1;
								}
							}
						}
						else if (restart && !cur_nod->gap_parent->active && reactive && (rew_pos == 0 || !cur_nod->gap_parent->del_psa)) {
							cur_nod->gap_parent->active = 1;
							cur_nod->gap_parent->del_psa = 0;
							++num_act_arcs;
						}
					}


					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					if (!cur_nod->children.empty()) {
						for (int j = cur_nod->children.size() - 1; j >= 0; --j) {
							if (init && (restart || cur_nod->children[j]->active)) {
								if (cur_nod->children[j]->child->parents.empty() || cur_nod->children[j]->reward[rew_pos] + cur_nod->children[j]->child->max_fr[rew_pos] + cur_nod->max_to[rew_pos] + 0.00001 < incumb - bound) {
									delete cur_nod->children[j];
									cur_nod->children.erase(cur_nod->children.begin() + j);
								}
							}
						}
						cur_nod->children.shrink_to_fit();
						cur_nod->parents.shrink_to_fit();
					}

					if ((num_inact_arc > 0 && num_inact_arc == cur_nod->parents.size() && (cur_nod->gap_parent == NULL || !cur_nod->gap_parent->active)) || (i > 0 && cur_nod->parents.empty() && (cur_nod->gap_parent == NULL || !cur_nod->gap_parent->active))) {
						cur_nod->active = 0;
						++num_inact_nod;
					}
					else
						cur_nod->active = 1;
				}

				num_act_nodes += (*it)->size() - num_inact_nod;

				if (num_inact_nod == (*it)->size())
					++num_inact_layer;

		}

			if (MDDs[s]->layers[i]->nodes.size() == num_inact_layer) {
				MDDs[s]->layers[i]->active = 0;
				cout << "Layer inactive\n";
				return 0;
			}
			else
				MDDs[s]->layers[i]->active = 1;
		}
	}

	cout << "Number of arcs after filtering: " << num_arcs << "   Number of active arcs: " << num_act_arcs << endl;

	if (!vars.empty() && ((!init && init_num > num_act_arcs) || (init && init_num < num_act_arcs + 1000000))) {
		cout << "deleting inactive vars\n";
		del_vars.endElements();
		vector<Var> new_vars;
		new_vars.reserve(num_act_arcs);
		int arc_num = 0;
		for (int i = 0; i < vars.size(); ++i) {
			if (!del[i]) {
				new_vars.emplace_back(vars[i]);
				vars[i].arc->ID = arc_num;
				++arc_num;
			}
		}
		vars.swap(new_vars);
	}

	del_vars.end();

	if (!init) {
		if (rew_pos == 0) {
			if (filt[0] == 0){
				filt[0] = num_act_nodes;
				filt[1] = num_act_arcs;
			}
			psa_feas = sol_exists(heur_sol);
			cout << "Heur solution exists: " << psa_feas << endl;
		}
		else {
			filt[4] = num_act_nodes;
			filt[5] = num_act_arcs;
			LP_feas = sol_exists(heur_sol);
			cout << "Heur solution exists: " << LP_feas << endl;
			if (LP_feas) {
				for (int s = strt; s < MDDs.size(); ++s) {
					for (int i = MDDs[s]->layers.size() - 1; i >= 0; --i) {
						for (vector<vector<Node*>*>::reverse_iterator it = MDDs[s]->layers[i]->nodes.rbegin(); it != MDDs[s]->layers[i]->nodes.rend(); ++it) {
							for (int n = (*it)->size() - 1; n >= 0; --n) {
								Node* cur_nod = (*it)->at(n);
								for (int j = 0; j < cur_nod->parents.size(); ++j) {
									if (!cur_nod->parents[j]->active)
										cur_nod->parents[j]->del_psa = 1;
								}
								if (cur_nod->gap_parent && !cur_nod->gap_parent->active)
									cur_nod->gap_parent->del_psa = 1;
							}
						}
					}
				}
			}
		}
	}
	else {
		psa_feas = 1;
		LP_feas = 1;
		if (rew_pos == 0 && filt[2] == 0) {
			filt[2] = num_act_nodes;
			filt[3] = num_act_arcs;
		}
		else {
			filt[6] = num_act_nodes;
			filt[7] = num_act_arcs;
		}
	}

	return 1;
}


void set_psa_del() {
	for (int s = 0; s < MDDs.size(); ++s) {
		for (int i = MDDs[s]->layers.size() - 1; i >= 0; --i) {
			for (vector<vector<Node*>*>::reverse_iterator it = MDDs[s]->layers[i]->nodes.rbegin(); it != MDDs[s]->layers[i]->nodes.rend(); ++it) {
				for (int n = (*it)->size() - 1; n >= 0; --n) {
					Node* cur_nod = (*it)->at(n);
					for (int j = 0; j < cur_nod->parents.size(); ++j) {
						if (!cur_nod->parents[j]->active && cur_nod->parents[j]->del_psa == 0){
							cur_nod->parents[j]->active = 1;
							++num_act_arcs;
						}
					}
					if (cur_nod->gap_parent && !cur_nod->gap_parent->active && cur_nod->gap_parent->del_psa == 0){
						cur_nod->gap_parent->active = 1;
						++num_act_arcs;
					}
				}
			}
		}
	}
}
