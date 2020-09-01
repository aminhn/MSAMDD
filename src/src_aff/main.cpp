#include <iostream>
#include <string.h>
#include <string>
#include "load_inst.hpp"
#include "build_mdd.hpp"
#include "benders.hpp"
//#include "utility.hpp"


using namespace std;

double opn_pen = 12, ext_pen = 2.22, incumb = -IloInfinity;
string item_file, subs_file;
vector<int> time_limit; 
vector<double> filt(8, 0);							//psa_filt_nod1, psa_filt_arc1, psa_filt_nod2, psa_filt_arc2, LP_filt_nod1, LP_filt_arc1, LP_filt_nod2, LP_filt_arc2
vector<double> elapsed_time(6, 0);						//overall, MDD, warmstart, heur, MIP, phase1
clock_t start_time;
IloEnv env;

int main(int argc, char* argv[]) {

	string output_file, start_file;

	int time_limit_in = 36000;

	subs_file = "./Data/blosum.ncbi";

	for (int i = 0; i < argc; ++i){
		if (argv[i][0] !='-' || isdigit(argv[i][1]))
			continue;
		else if (strcmp(argv[i], "-in") == 0) {
			item_file = argv[i + 1];
			size_t lastindex = item_file.find_last_of(".");
			if (lastindex > 0)
				item_file = item_file.substr(0, lastindex);
		}
		else if (strcmp(argv[i], "-out") == 0)
			output_file = argv[i + 1];
		else if (strcmp(argv[i], "-start") == 0)
			start_file = argv[i + 1];
		else if (strcmp(argv[i], "-time") == 0)
			time_limit_in = stoi(argv[i + 1]);
		else if (strcmp(argv[i], "-op") == 0)
			opn_pen = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-ep") == 0)
			ext_pen = stod(argv[i + 1]);
                else if (strcmp(argv[i], "-submat") == 0)
                        subs_file = argv[i + 1];
		else 
			cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!Command " << argv[i] << " not recognized and skipped.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	}

	cout << "\n\n-------------------------------------------------------------------------- " << item_file << "--------------------------------------------------------------------------\n\n";

	time_limit.push_back(time_limit_in);					//overall
	time_limit.push_back(time_limit_in);					//warmstart
	time_limit.push_back(3600);						//LP limit


	start_time = clock();

	if (!Load_subsmat()) {
		cout << "File invalid, exiting.\n";
		cin.get();
		return 1;
	}

	if (!Load_instance(item_file)) {
		cout << "File invalid, exiting.\n";
		cin.get();
		return 1;
	}

	if (!start_file.empty()) {
		cout << "Loading starting solution:\n";
		incumb = Load_sol(start_file);
		cout << "Starting solution gave an incumbent of " << incumb << endl;
	}

	cout << "Instance loaded\n";

	clock_t kk = clock();

	double incumb_muscle = run_MUSCLE(item_file, 0);

	if (incumb_muscle > incumb)
		incumb = incumb_muscle;

	elapsed_time[3] += give_time(clock() - kk);

	cout << "heuristic solution value: " << incumb_muscle << endl;

	kk = clock();

	Build_mdd();

	elapsed_time[1] += give_time(clock() - kk);

	cout << "MDDs built\n";

	IloCplex cplex(env);

	kk = clock();

	init_warm(0, cplex);

	elapsed_time[2] += give_time(clock() - kk);

	Benders_IP(cplex);

	elapsed_time[5] = elapsed_time[0];

	if (time_limit[0] > elapsed_time[0])
		cout << "Benders algorithm converged\n\n";
	else
		cout << "Phase 1 time limit reached\n\n";

	double ph1_val = incumb;

	cout << "Phase 1 total time: " << elapsed_time[0] << endl;

	if (elapsed_time[0] < time_limit[0]) {
		
		kk = clock();

		init_warm(2, cplex);

		elapsed_time[2] += give_time(clock() - kk);

		Benders_IP(cplex);

		elapsed_time[0] = give_time(clock() - start_time);

		cout << "Phase 2 total time: " << give_time(clock() - start_time) << endl;

	}


	env.end();

	cout << "\nEND\n\n";

	string del_file = item_file + "_s.txt";
	remove(del_file.c_str());

	cout << "Total time: " << elapsed_time[0] << "  Phase 1: " << elapsed_time[5] << "  MDD: " << elapsed_time[1] << "  Warmstart: " << elapsed_time[2] << "  Heur: " << elapsed_time[3] << "  MIP: " << elapsed_time[4] << endl;
	cout << "Gap is: " << gap << endl;

	ofstream file3;
	file3.open(output_file);
	for (int i = 0; i < N; ++i) {
		file3 << seq_names[i] << "\n";
		for (vector<int>::iterator it = heur_sol[i]->begin(); it != heur_sol[i]->end(); ++it) {
			file3 << (char)(*it);
		}
		file3 << endl;
	}

	file3.close();

}