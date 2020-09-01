#include "load_inst.hpp"
#include "utility.hpp"
#include <map>

vector<vector<int>*> strings, heur_sol, subsmat;
vector<int>  nod_to_str, alphabet;
vector<int> comul_siz, comul_mdd;
map<int, int> subsref;

map<string, int> seq_name_ord;

vector<string> seq_names;

size_t N = 0, M = 0, tot_chars = 0;
int num_arcs = 0, num_act_arcs = 0;

bool Load_instance(string &inst) {

	ifstream file(inst + ".fa");

	if (file.good()) {

		string line, seq;
		int lin_num = 0, seq_num = 0;

		while (getline(file, line)) {
			if (line.empty())
				continue;
			if (line[0] != '>') {
				seq += line;
				++lin_num;
			}
			else {
				seq_names.emplace_back(line);
				if (lin_num != 0) {
					strings.push_back(new vector<int>);
					for (string::iterator it = seq.begin(); it != seq.end(); it++)
						strings.back()->push_back((int)*it);
					seq.clear();
				}
			}
		}
		strings.push_back(new vector<int>);
		for (string::iterator it = seq.begin(); it != seq.end(); it++)
			strings.back()->push_back((int)*it);
	}
	else {
		cout << "!!!!!! No such file exists: " << inst << " !!!!!!\n";
		return 0;
	}

	N = strings.size();

	vector<size_t> idx(N);
	vector<int> temp_siz(N);
	for (int i = 0; i < N; ++i) {
		temp_siz[i] = strings[i]->size();
		idx[i] = i;
	}

	sort(idx.begin(), idx.end(), [&temp_siz](size_t i, size_t j) {return temp_siz[i] > temp_siz[j]; });

	vector<vector<int>*> temp_strings(N);
	vector<string> temp_names(N);

	for (int i = 0; i < N; ++i) {
		temp_strings[i] = strings[idx[i]];
		temp_names[i] = seq_names[idx[i]];
		seq_name_ord[temp_names[i]] = i;
	}

	strings.swap(temp_strings);
	seq_names.swap(temp_names);	

	for (int i = 0; i < strings.size(); i++) {
		if (strings[i]->size() > M)
			M = strings[i]->size();
	}

	comul_siz.reserve(N + 1);

	size_t mdd_count = 0;
	for (int i = 0; i < strings.size(); ++i) {
		comul_mdd.push_back(mdd_count);
		mdd_count += strings.size() - i - 1;
		comul_siz.push_back(tot_chars);
		tot_chars += strings[i]->size();
	}
	comul_siz.push_back(tot_chars);

	nod_to_str.reserve(comul_siz.back());
	for (int i = 0; i < strings.size(); ++i) {
		for (int j = 0; j < strings[i]->size(); ++j) {
			nod_to_str.push_back(i);
		}
	}

	return 1;

}

double Load_sol(string& inst) {
	
	ifstream file2(inst);
	vector<int> ord(N);
	vector<vector<int>*> temp(N);
	string name;

	if (file2.good()) {

		string line, seq;
		int lin_num = 0;

		while (getline(file2, line)) {
			if (line.empty())
				continue;
			if (line[0] != '>') {
				seq += line;
				++lin_num;
			}
			else{
				if (lin_num != 0) {
					size_t pl = seq_name_ord[name];
					temp[pl] = new vector<int>;
					for (string::iterator it = seq.begin(); it != seq.end(); it++)
						temp[pl]->emplace_back((int)*it);
					seq.clear();
				}
				name = line;
			}
		}
		size_t pl = seq_name_ord[name];
		temp[pl] = new vector<int>;
		for (string::iterator it = seq.begin(); it != seq.end(); it++)
			temp[pl]->push_back((int)*it);
	}
	else {
		cout << "!!!!!! No solution file exists: " << inst << " !!!!!!\n";
		return -IloInfinity;
	}

	double val = get_value(temp);

	if (heur_sol.size() == 0 || val > incumb) {
		disp_sol(temp);
		for (vector<vector<int>*>::iterator it = heur_sol.begin(); it != heur_sol.end(); ++it)
			delete (*it);
		heur_sol.swap(temp);
		check_feas(heur_sol);
	}
	else
		for (vector<vector<int>*>::iterator it = temp.begin(); it != temp.end(); ++it)
			delete (*it);

	return val;

}


bool Load_subsmat() {

	ifstream file(subs_file);

	if (file.good()) {

		string line, seq;
		
		int lin_num = 0, col_num, ditem;

		while (getline(file, line)){
			col_num = 0;
			istringstream word(line);
			string itm;
			if (lin_num > 0)
				subsmat.push_back(new vector<int>);
			while (word >> itm) {
				if (lin_num == 0) {
					subsref[(int)itm[0]] = col_num;
					alphabet.push_back((int)itm[0]);
				}
				else if (col_num > 0) {
					ditem = stoi(itm);
					subsmat.back()->push_back(ditem);
				}
				++col_num;
			}
			++lin_num;
		}

	}
	else {
		cout << "!!!!!! Substitution matrix file: " << subs_file << " does not exist!!!!!!\n";
		return 0;
	}

	return 1;

}


void write_sol(list<vector<int>>& matrix, string& inst) {

	ofstream file;
	string wrt = inst + ".fa";
	file.open(wrt);

	for (int i = 0; i < N; ++i) {
		file << seq_names[i] << "\n";
		for (list<vector<int>>::iterator it = matrix.begin(); it != matrix.end(); ++it) {
			if (it->at(i) == -1)
				file << "-";
			else
				file << (char)strings[i]->at(it->at(i));
		}
		file << endl;
	}

	file.close();

}