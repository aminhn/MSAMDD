#include "heur.hpp"

void make_feas(list<vector<int>>& feas_sol);
void fix_sol(list<vector<int>>& feas_sol, vector<vector<bool>>& placed);

void heuristic() {

	list<vector<int>> feas_sol;

	make_feas(feas_sol);

	string heur_f = item_file + "_h";

	write_sol(feas_sol, heur_f);

	double val = run_MUSCLE(heur_f, 1);

	cout << "Heuristic solution value is: " << val << endl;

	if (val > incumb) {
		cout << "Incumbent found in heuristic!\n";
		incumb = val;
	}
}


void find_incumb(IloCplex& cplex) {

	int num_sol = cplex.getNMIPStarts();

	cout << "Looking for feasible solutions in " << num_sol << " MIP starts\n";

	IloModel model = cplex.getModel();

	IloNumVarArray mip_vars(env);
	for (int i = 0; i < vars.size(); ++i)
		mip_vars.add(vars[i].var);

	for (int i = num_sol - 2; i >= 0; --i) {
		IloNumArray mip_vals(env);
		IloBoolArray mip_bin(env);
		cplex.getMIPStart(i, mip_vars, mip_vals, mip_bin);
		for (int j = 0; j < vars.size(); ++j) {
			vars[j].val = mip_vals[j];
		}

		IloRange unb_ray;
		if (!build_SP(unb_ray)) {
			model.add(unb_ray);
			double prev_incumb = incumb;
			heuristic();
		}
		else {
			unb_ray.end();
			cout << "A FEASIBLE SOLUTION FOUND FROM MIP STARTS!\n";
		}
		mip_vals.end();
		mip_bin.end();
	}

	mip_vars.end();
}


void make_feas(list<vector<int>>& feas_sol) {

	vector<double> cost;
	vector<int> act_sol;
	vector<size_t> idx;

	cost.reserve(vars.size());
	act_sol.reserve(vars.size());
	idx.reserve(vars.size());

	int num_sol = 0;
	for (int i = 0; i < vars.size(); ++i) {
		if (vars[i].val > 0 && vars[i].arc->parent->state != vars[i].arc->child->state) {
			cost.push_back(vars[i].arc->reward[0]);
			act_sol.push_back(i);
			idx.push_back(num_sol);
			++num_sol;
		}
	}

	sort(idx.begin(), idx.end(), [&cost](size_t i, size_t j) {return cost[i] > cost[j]; });

	vector<vector<bool>> placed;
	for (int i = 0; i < N; ++i)
		placed.emplace_back(vector<bool>(strings[i]->size(), 0));

	int stri = MDDs[vars[act_sol[idx[0]]].arc->parent->mdd_count]->seq_fr, strj = MDDs[vars[act_sol[idx[0]]].arc->parent->mdd_count]->seq_to;
	int let1 = strings[stri]->size() - vars[act_sol[idx[0]]].arc->parent->lvl - 1;
	int let2 = vars[act_sol[idx[0]]].arc->child->state - 1;
	vector<int> new_col(N, -1);
	new_col[stri] = let1;
	new_col[strj] = let2;
	feas_sol.emplace_back(new_col);
	placed[stri][let1] = 1;
	placed[strj][let2] = 1;

	for (int i = 1; i < num_sol; ++i) {
		int stri = MDDs[vars[act_sol[idx[i]]].arc->parent->mdd_count]->seq_fr, strj = MDDs[vars[act_sol[idx[i]]].arc->parent->mdd_count]->seq_to;
		int let1 = strings[stri]->size() - vars[act_sol[idx[i]]].arc->parent->lvl - 1;
		int let2 = vars[act_sol[idx[i]]].arc->child->state - 1;

		if (placed[stri][let1] && placed[strj][let2])
			continue;

		//cout << "##### " << i << " " << stri + 1 << " " << strj + 1 << " " << let1 << " " << let2 << endl;

		int pos = 0;
		for (list<vector<int>>::iterator it = feas_sol.begin(); it != feas_sol.end(); ++it) {
			if (it->at(stri) < let1 && it->at(strj) < let2) {
				++pos;
				if (next(it, 1) == feas_sol.end()) {
					//	cout << "1\n";
					vector<int> new_col(N, -1);
					new_col[stri] = let1;
					new_col[strj] = let2;
					feas_sol.emplace_back(new_col);
					placed[stri][let1] = 1;
					placed[strj][let2] = 1;
				}
			}
			else if (it->at(stri) == let1) {
				//cout << "2\n";
				bool valid = 1;
				for (int j = 1; next(it, j) != feas_sol.end(); ++j) {
					if (next(it, j)->at(strj) != -1) {
						if (let2 >= next(it, j)->at(strj))
							valid = 0;
						break;
					}
				}
				if (valid && it->at(strj) == -1) {
					it->at(strj) = let2;
					placed[strj][let2] = 1;
				}
				break;
			}
			else if (it->at(strj) == let2) {
				//cout << "3\n";
				bool valid = 1;
				for (int j = 1; next(it, j) != feas_sol.end(); ++j) {
					if (next(it, j)->at(stri) != -1) {
						if (let1 >= next(it, j)->at(stri))
							valid = 0;
						break;
					}
				}
				if (valid && it->at(stri) == -1) {
					it->at(stri) = let1;
					placed[stri][let1] = 1;
				}
				break;
			}
			else if (it->at(stri) > let1) {
				//cout << "3\n";
				bool valid = 1;
				for (int j = 0; next(it, j) != feas_sol.end(); ++j) {
					if (next(it, j)->at(strj) != -1) {
						if (let2 >= next(it, j)->at(strj))
							valid = 0;
						break;
					}
				}
				if (valid) {
					vector<int> new_col(N, -1);
					new_col[stri] = let1;
					new_col[strj] = let2;
					feas_sol.insert(it, new_col);
					placed[stri][let1] = 1;
					placed[strj][let2] = 1;
				}
				break;
			}
			else if (it->at(strj) > let2) {
				//cout << "4\n";
				bool valid = 1;
				for (int j = 0; next(it, j) != feas_sol.end(); ++j) {
					if (next(it, j)->at(stri) != -1) {
						if (let1 >= next(it, j)->at(stri))
							valid = 0;
						break;
					}
				}
				if (valid) {
					vector<int> new_col(N, -1);
					new_col[stri] = let1;
					new_col[strj] = let2;
					feas_sol.insert(it, new_col);
					placed[stri][let1] = 1;
					placed[strj][let2] = 1;
				}
				break;
			}
		}
		/*for (int i = 0; i < N; ++i) {
		for (list<vector<int>>::iterator it = feas_sol.begin(); it != feas_sol.end(); ++it) {
		if (it->at(i) != 45 && it->at(i) / 10 < 1)
		cout << it->at(i) << "  ";
		else
		cout << it->at(i) << " ";
		}
		cout << endl;
		}
		cin.get();*/
	}

	fix_sol(feas_sol, placed);

}

void fix_sol(list<vector<int>>& feas_sol, vector<vector<bool>>& placed) {

	cout << "Fixing heuristic solution\n";

	for (int i = 0; i < N; ++i) {
		int pos = 0;
		for (list<vector<int>>::iterator it = feas_sol.begin(); it != feas_sol.end(); ++it) {
			if (it->at(i) == -1)
				continue;
			else if (!placed[i][pos]) {
				int last_g = 1;
				if (it != feas_sol.begin()) {
					while (prev(it, last_g)->at(i) == -1)
						++last_g;
				}
				if (last_g == 1) {
					vector<int> new_col(N, -1);
					new_col[i] = pos;
					feas_sol.insert(it, new_col);
					placed[i][pos] = 1;
					--it;
				}
				else {
					prev(it, last_g - 1)->at(i) = pos;
					placed[i][pos] = 1;
					--it;
				}
			}
			++pos;
		}
		if (pos != placed[i].size()) {
			while (pos != placed[i].size()) {
				vector<int> new_col(N, -1);
				new_col[i] = pos;
				feas_sol.emplace_back(new_col);
				placed[i][pos] = 1;
				++pos;
			}
		}
	}
	cout << "Fin\n";
}




