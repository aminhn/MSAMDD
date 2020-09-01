#include "benders.hpp"

list<IloRange> viol_cons;
double gap, LB_guess_psa, LB_guess_LP, stp_siz, UB;
bool psa_feas = 0, LP_feas = 0;

void init_warmstart(int update) {

	clock_t kk = clock();

	UB = IloInfinity;

	if (update == 0) {

		stp_siz = (opt_psa - incumb) / ((double)N * 75 / 4 - 50);

		cout << "step size: " << stp_siz << endl;

		cout << "Lngst path time: " << give_time(clock() - kk) << endl;
		cout << "Optimal PSA " << opt_psa << endl;

		LB_guess_psa = opt_psa - stp_siz, LB_guess_LP = LB_guess_psa;
		//LB_guess_psa = incumb, LB_guess_LP = LB_guess_psa;
		cout << "LB_guess " << LB_guess_LP << endl;
	}
	else if (update == 1) {
		if (!psa_feas){
			LB_guess_psa -= stp_siz;
			filt[0] = 0;
		}
		if (!LP_feas)
			LB_guess_LP -= stp_siz;
		cout << "new guess: " << LB_guess_psa << " " << LB_guess_LP << endl;
	}
	else {
		LB_guess_psa = incumb, LB_guess_LP = LB_guess_psa;
		filt[2] = 0;
		psa_feas = 1;
		LP_feas = 1;
		cout << "Guaranteed LB " << LB_guess_LP << endl;
	}


}

void define_model(IloCplex& cplex, bool restart) {
	cout << "Building cplex model\n";

	IloModel model(env);

	build_model(model);

	cplex = IloCplex(model);

	cout << "Cplex model Built\n";

	opt_psa = 0;
	for (vector<MDD*>::iterator it = MDDs.begin(); it != MDDs.end(); ++it) {
		(*it)->comp_lng_to(0, restart);
		opt_psa += (*it)->lng_path->max_to[0];
		(*it)->comp_lng_fr(0, restart);
	}
	cout << endl;

	cout << "NEW calculated Optimal PSA " << opt_psa << endl;

	for (vector<MDD*>::iterator it = MDDs.begin(); it != MDDs.end(); ++it)
		(*it)->extract_lngp();

	viol_cons.clear();
}

void init_warm(int update, IloCplex& cplex) {

	cout << "warmstarting\n\n";

	clock_t kk = clock();

	init_warmstart(update);

	bool infeas = 1, first_iter = 1;
	while (infeas && (give_time(clock() - start_time) < time_limit[1] || (update == 2 && give_time(clock() - start_time) < time_limit[0]))) {

		cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INITIALIZE\n";

		env.end();
		env = IloEnv();
		vars.clear();

		bool restart = !(update == 1 && !first_iter && psa_feas);

		for (vector<MDD*>::iterator it = MDDs.begin(); it != MDDs.end(); ++it) {
			(*it)->comp_lng_to(0, restart);
			(*it)->comp_lng_fr(0, restart);
		}
		cout << endl;

		if (!restart) {
			filter(opt_psa, LB_guess_psa, 0, 0, 1, 0);
		}
		else {
			while (!filter(opt_psa, LB_guess_psa, 0, 0, 1, 1) || (update == 1 && !psa_feas)) {
				LB_guess_psa -= stp_siz;
				LB_guess_LP = LB_guess_psa;
				filt[0] = 0;
				cout << "PSA guess invalid, new guess: " << LB_guess_psa << endl;
			}
			first_iter = 0;
		}


		define_model(cplex, restart);

		if (update == 2) {
			if (num_act_arcs > 30000000){
				cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
				cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
			}
			cplex.setParam(cplex.RootAlg, 3);
			cplex.setParam(cplex.Threads, 1);
			//cplex.setParam(IloCplex::Param::TimeLimit, time_limit[2]);
			cout << "Solving LP\n";
			cplex.solve();
			cplex.setParam(cplex.RootAlg, 0);
		}


		IloModel model = cplex.getModel();

		double prev_val = opt_psa + 10, LP_val;
		UB = IloInfinity;
		bool UB_updt = 1;
		while (give_time(clock() - start_time) < time_limit[1] || (update == 2 && give_time(clock() - start_time) < time_limit[0])) {

			elapsed_time[0] = give_time(clock() - start_time);

			cout << "\nElapsed time: " << elapsed_time[0] << endl << endl;

			if (UB_updt) {
				cout << "GETTING VIOL\n";
				list<IloRange> cur_viol;
				get_viol(cur_viol);
	
				cout << "NUM VIO " << cur_viol.size() << endl;

				for (list<IloRange>::iterator it = cur_viol.begin(); it != cur_viol.end(); ++it) {
					model.add(*it);
					viol_cons.emplace_back(*it);
				}
	
				cout << "Overall violated size: " << viol_cons.size() << endl;
			}

			cplex.setParam(cplex.Threads, 1);
			cplex.setParam(IloCplex::Param::TimeLimit, time_limit[2]);
			//cplex.setOut(env.getNullStream());
			if (num_act_arcs > 30000000){
				cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
				cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
			}
			else{
				cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 1);
				cplex.setParam(IloCplex::Param::Emphasis::Memory, 0);
			}

			kk = clock();
			cout << "Solving LP\n";
			cplex.solve();
			cout << "Solve Time: " << give_time(clock() - kk) << endl;

			if (cplex.getStatus() == IloAlgorithm::Infeasible || cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded) {
				cout << "cplex infeasible\n";
				if (!psa_feas) {
					LB_guess_psa -= stp_siz;
					LB_guess_LP = LB_guess_psa;
					filt[0] = 0;
					cout << "PSA guess invalid, new guess: " << LB_guess_LP << endl;
					break;
				}
				LB_guess_LP -= stp_siz;
				cout << "LB guess invalid, new guess: " << LB_guess_LP << endl;
				break;
			}

			LP_val = cplex.getObjValue();
			UB_updt = 0;

			if (LP_val < UB && cplex.getCplexStatus() != IloCplex::AbortTimeLim) {
				UB = LP_val;
				UB_updt = 1;
			}

			if (UB < LB_guess_LP && update != 2) {
				LB_guess_LP = UB - stp_siz;
				cout << "LP guess readjusted to optimal objective - stp\n";
			}

			cout << "stat: " << cplex.getStatus() << " " << cplex.getCplexStatus() << " obj: " << cplex.getObjValue() << endl;

			for (size_t i = 0; i < vars.size(); ++i) {
				vars[i].val = round(cplex.getValue(vars[i].var) * 10000) / 10000;
				if (vars[i].val > 0) {
					vars[i].arc->reward[1] = 0;
				}
				else
					vars[i].arc->reward[1] = round(cplex.getReducedCost(vars[i].var) * 10000) / 10000;
			}

			if (prev_val - UB < 0.1 && (UB_updt || cplex.getStatus() == IloAlgorithm::Optimal)) {

				cout << "\n\n num viol before: " << viol_cons.size() << endl;

				IloRangeArray del_cons(env);
				list<IloRange>::iterator it = viol_cons.begin();
				while (!viol_cons.empty() && it != viol_cons.end()) {
					if (cplex.getDual(*it) == 0) {
						del_cons.add(*it);
						it = viol_cons.erase(it);
					}
					else{
						++it;
					}
				}

				cout << "deleting...\n";

				del_cons.endElements();
				del_cons.end();

				cout << "\n\nOverall violated size: " << viol_cons.size() << endl << endl;
			}

			for (vector<MDD*>::iterator it = MDDs.begin(); it != MDDs.end(); ++it) {
				(*it)->comp_lng_to(1, 0);
				(*it)->comp_lng_fr(1, 0);
			}
			cout << endl;

			int prev_num_act = num_act_arcs;

			if (!filter(opt_psa, LB_guess_LP, LP_val, 1, 0, 0)) {
				cout << "LP guess invalid\n";
				if (!psa_feas) {
					LB_guess_psa -= stp_siz;
					LB_guess_LP = LB_guess_psa;
					filt[0] = 0;
				}
				else {
					LB_guess_LP -= stp_siz;
					set_psa_del();
				}
				break;
			}

			if (cplex.getCplexStatus() == IloCplex::AbortTimeLim && prev_num_act == num_act_arcs) {
				time_limit[2] += 1800;
				cout << "No Improvement in bound, readjusting LP time limit to: " << time_limit[2] << endl;
				continue;
			}
			else if (cplex.getCplexStatus() == IloCplex::AbortTimeLim && UB < LP_val) {
				UB = LP_val;
				UB_updt = 1;
			}


			if (update == 2 && prev_num_act > 1000000 + num_act_arcs)
				break;

			if (update == 1 && !LP_feas) {
				cout << "LP guess invalid in update\n";
				LB_guess_LP -= stp_siz;
				cout << "new guess " << LB_guess_LP << endl;
				set_psa_del();
				break;
			}

			if (prev_val - UB < 0.1 && cplex.getCplexStatus() != IloCplex::AbortTimeLim) {
				cout << "No Improvement in bound " << UB << " " << prev_val - UB << endl;
				infeas = 0;
				cout << "GETTING VIOL\n";
				list<IloRange> cur_viol;
				get_viol(cur_viol);

				cout << "NUM VIO " << cur_viol.size() << endl;

				for (list<IloRange>::iterator it = cur_viol.begin(); it != cur_viol.end(); ++it) {
					model.add(*it);
					viol_cons.emplace_back(*it);
				}

				cout << "Overall violated size: " << viol_cons.size() << endl;

				break;	
			}

			double prev_incumb = incumb;
			IloRange unb_ray;
			if (!build_SP(unb_ray))
				model.add(unb_ray);
			else {
				unb_ray.end();
				cout << "FEASIBILE SOLUTION!!\n";
				if (incumb != prev_incumb) {
					if (update == 2)
						LB_guess_psa = incumb, LB_guess_LP = LB_guess_psa;
					if (update != 0)
						break;
				}
			}

			if (UB_updt)
				prev_val = UB;

		}
	}

	elapsed_time[0] = give_time(clock() - start_time);

	cout << "\nElapsed time: " << elapsed_time[0] << endl << endl;

	gap = (UB - incumb) / abs(incumb) * 100;

	cout << "Gap is: " << gap << endl;

	cout << "\n\n---------------------------------FILTER DONE-----------------------------------\n\n";

}


void Benders_IP(IloCplex& master) {

	cout << "\n\nStarting Benders Decomposition\n\n";

	elapsed_time[0] = give_time(clock() - start_time);

	IloModel model = master.getModel();

	for (vector<Var>::iterator it = vars.begin(); it != vars.end(); ++it)
		model.add(IloConversion(env, it->var, ILOINT));

	clock_t kk;

	gap = (UB - incumb) / abs(incumb) * 100;

	if (gap < 0.01) {
		cout << "Upper bound at incumbent, gap = 0, UB = " << UB << " Inc = " << incumb << endl;
		gap = 0;
		disp_sol(heur_sol);
	}

	IloNumVarArray mipvars(env);
	IloNumArray mipvals(env);
	bool add_mipstr = 0;
	if (psa_feas && LP_feas) {
		cout << "Adding MIP start from heuristic solution\n";
		if (!get_MIPstart(mipvars, mipvals))
			cout << "HEURISTIC not in MDDs\n";
		else
			add_mipstr = 1;
	}

	double IP_val;

	while (gap >= 0.01 && elapsed_time[0] < time_limit[0]) {

		cout << "\nElapsed time: " << elapsed_time[0] << endl << endl;

		if (add_mipstr)
			master.addMIPStart(mipvars, mipvals, IloCplex::MIPStartNoCheck);
		else
			cout << "No MIP start added\n";

		kk = clock();
		master.setParam(master.Threads, 1);
		master.setParam(IloCplex::Param::TimeLimit, 100 + time_limit[0] - (give_time(clock() - start_time)));
		//master.setOut(env.getNullStream());
		cout << "Solving MP\n";
		master.solve();
		elapsed_time[4] += give_time(clock() - kk);

		if (master.getStatus() == IloAlgorithm::Infeasible || master.getStatus() == IloAlgorithm::InfeasibleOrUnbounded) {
			cout << "Benders master infeasible\n";
			kk = clock();
			init_warm(1, master);
			elapsed_time[2] += give_time(clock() - kk);
			model = master.getModel();
			for (vector<Var>::iterator it = vars.begin(); it != vars.end(); ++it)
				model.add(IloConversion(env, it->var, ILOINT));
			mipvars = IloNumVarArray(env);
			mipvals = IloNumArray(env);
			if (!get_MIPstart(mipvars, mipvals)) {
				cout << "ERROR: HEURISTIC not in MDDs, even after making sure!\n";
				cin.get();
				break;
			}
			else
				add_mipstr = 1;
			continue;
		}
		else if (master.getStatus() != IloAlgorithm::Optimal)
			cout << "TIME LIMIT REACHED in CPLEX\n";


		IP_val = master.getBestObjValue();

		if (IP_val < UB)
			UB = IP_val;

		gap = (UB - incumb) / abs(incumb) * 100;

		cout << "Gap is: " << gap << endl;
		cout << "STAT: " << master.getStatus() << " Solve Time: " << give_time(clock() - kk) << " obj: " << IP_val << endl;

		if (gap < 0.01) {
			cout << "Gap reached, optimal solution found!\n";
			disp_sol(heur_sol);
			elapsed_time[0] = give_time(clock() - start_time);
			break;
		}

		if (master.getStatus() != IloAlgorithm::Optimal && (master.getStatus() != IloAlgorithm::Feasible || master.getNMIPStarts() == 0)) {
			cout << "No solution found by Cplex\n";
			elapsed_time[0] = give_time(clock() - start_time);
			break;
		}


		for (size_t i = 0; i < vars.size(); ++i) {
			vars[i].val = round(master.getValue(vars[i].var));
		}

		cout << "GETTING VIOL\n";
		list<IloRange> cur_viol;
		get_viol(cur_viol);

		cout << "NUM VIO " << cur_viol.size() << endl;

		for (list<IloRange>::iterator it = cur_viol.begin(); it != cur_viol.end(); ++it) {
			model.add(*it);
			viol_cons.emplace_back(*it);
		}

		cout << "Overall violated size: " << viol_cons.size() << endl;

		IloRange unb_ray;
		if (!build_SP(unb_ray)) {
			model.add(unb_ray);
			double prev_incumb = incumb;
			kk = clock();
			heuristic();
			elapsed_time[3] += give_time(clock() - kk);
			if (prev_incumb != incumb && psa_feas && LP_feas) {
				gap = (UB - incumb) / abs(incumb) * 100;
				cout << "Gap is: " << gap << endl;
				mipvars.end();
				mipvals.end();
				mipvars = IloNumVarArray(env);
				mipvals = IloNumArray(env);
				if (!get_MIPstart(mipvars, mipvals) && incumb >= UB) {
					cout << "New heuristic not in MDDs, ending phase 1\n";
					elapsed_time[0] = give_time(clock() - start_time);
					break;
				}
				else
					add_mipstr = 1;
			}
			else if (prev_incumb != incumb){
				gap = (UB - incumb) / abs(incumb) * 100;
				cout << "Gap is: " << gap << endl;
			}
		}
		else if ((UB - incumb) / abs(incumb) * 100 < 0.01){
			unb_ray.end();
			cout << "OPTIMAL SOLUTION FOUND!\n";
			elapsed_time[0] = give_time(clock() - start_time);
			gap = 0;
			break;
		}
		else
			unb_ray.end();

		double prev_incumb = incumb;

		kk = clock();
		find_incumb(master);
		elapsed_time[3] += give_time(clock() - kk);

		if (prev_incumb != incumb && LP_feas && psa_feas) {
			gap = (UB - incumb) / abs(incumb) * 100;
			cout << "Gap is: " << gap << endl;
			mipvars.end();
			mipvals.end();
			mipvars = IloNumVarArray(env);
			mipvals = IloNumArray(env);
			if (!get_MIPstart(mipvars, mipvals) && incumb >= UB) {
				cout << "New heuristic not in MDDs, ending phase 1\n";
				elapsed_time[0] = give_time(clock() - start_time);
				break;
			}
			else
				add_mipstr = 1;
		}
		else if (prev_incumb != incumb){
			gap = (UB - incumb) / abs(incumb) * 100;
			cout << "Gap is: " << gap << endl;
		}

		elapsed_time[0] = give_time(clock() - start_time);

	}

	mipvars.end();
	mipvals.end();

}