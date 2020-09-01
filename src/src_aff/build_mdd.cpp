#include "build_mdd.hpp"
#include "filter.hpp"
#include "build_model.hpp"

vector<MDD*> MDDs;
vector<vector<int> > value_ref;
vector<vector<double> > ext_gap_ref;
vector<vector<vector<psa_nod> >*> psa_mats;

double opt_psa = 0;

void get_psa(vector<int> seq1, vector<int> seq2);

void Build_mdd() {

	get_ext_ref(ext_gap_ref);

	/*for (int i = 0; i < ext_gap_ref.size(); ++i) {
		for (int j = 0; j < ext_gap_ref[i].size(); ++j) {
			cout << ext_gap_ref[i][j] << " ";
		}
		cout << endl;
	}

	cin.get();*/

	MDDs.reserve(strings.size() * (strings.size() - 1) / 2);
	psa_mats.reserve(MDDs.size());


	for (size_t i = 0; i < strings.size() - 1; ++i) {
		for (size_t j = i + 1; j < strings.size(); ++j) {
			value_ref = vector<vector<int>>(alphabet.size());
			get_psa(*strings[i], *strings[j]);
		}
	}

	cout << "Optimal PSA " << opt_psa << endl;

	int mdd_count = 0;
	for (size_t i = 0; i < strings.size() - 1; ++i) {
		for (size_t j = i + 1; j < strings.size(); ++j) {
			clock_t kk = clock();
			MDDs.emplace_back(new MDD(strings[i]->size(), i, j));
			cout << "S1: " << strings[i]->size() << " S2: " << strings[j]->size() << endl;
			value_ref = vector<vector<int>>(alphabet.size());

			//get_psa(*strings[i], *strings[j]);

			MDDs.back()->create_layer(strings[j]->size(), 0, mdd_count, 0);											//create root layer
			MDDs.back()->create_layer(strings[j]->size(), 1, mdd_count, 0);											//create lvl = 1 layer

			int pos = subsref[strings[i]->back()];
			get_align_values(*strings[j], value_ref[pos], pos);	
																												//create first layer arcs
			MDDs.back()->create_arcs_r(MDDs.back()->layers[0]->nodes, MDDs.back()->layers[1]->nodes, pos, MDDs.size() - 1);

			for (size_t s = 2; s < strings[i]->size(); ++s) {													//create layers 2 -> |S1| - 1 layers and arcs
			
				//cout << "layer " << s << " Built\n";
				
				pos = subsref[strings[i]->at(strings[i]->size() - s)];
				if (value_ref[pos].empty()) 
					get_align_values(*strings[j], value_ref[pos], pos);
				
				MDDs.back()->create_layer(strings[j]->size(), s, mdd_count, 0);
				MDDs.back()->create_arcs(MDDs.back()->layers[s - 1]->nodes, MDDs.back()->layers[s]->nodes, s - 1, pos, MDDs.size() - 1);

				//cout << "arcs built\n";
			}

			//cout << "last layer " << strings[i]->size() << " Built\n";

			pos = subsref[strings[i]->at(0)];
			if (value_ref[pos].empty())
				get_align_values(*strings[j], value_ref[pos], pos);

			MDDs.back()->create_layer(strings[j]->size(), strings[i]->size(), mdd_count, 1);
			MDDs.back()->create_arcs_t(MDDs.back()->layers[strings[i]->size() - 1]->nodes, MDDs.back()->layers.back()->nodes, strings[i]->size() - 1, pos, MDDs.size() - 1);
			//cout << "last arcs built\n";

			cout << give_time(clock() - kk) << endl;

			cout << "....................MDD...........................\n";

			MDDs.back()->comp_lng_fr(0, 1);

			num_act_arcs = num_arcs;

			++mdd_count;

			delete psa_mats[MDDs.size() - 1];

			//cin.get();

		}
	}

	num_act_arcs = num_arcs;

}

void get_psa(vector<int> seq2, vector<int> seq1) {

	psa_mats.emplace_back(new vector<vector<psa_nod> >(seq1.size() + 1));
	psa_mats.back()->at(0).reserve(seq2.size() + 1);
	for (size_t i = 0; i < seq2.size(); ++i) {
		psa_mats.back()->at(0).emplace_back(ext_gap_ref[0][i], i, ext_gap_ref[1][i], i);
	}
	psa_mats.back()->at(0).emplace_back(ext_gap_ref[0][seq2.size()], seq2.size(), ext_gap_ref[0][seq2.size()], seq2.size());

	for (size_t i = 1; i < psa_mats.back()->size(); ++i) {
		psa_mats.back()->at(i).reserve(seq2.size() + 1);
		psa_mats.back()->at(i).emplace_back(ext_gap_ref[0][i], i, ext_gap_ref[1][i], i);
	}

	double lr, ud, pes_val;
	for (size_t i = 1; i < seq1.size() + 1; ++i) {
		int pos = subsref[seq1[i - 1]];
		if (value_ref[pos].empty())
			get_align_values(seq2, value_ref[pos], pos);

		for (size_t j = 1; j < seq2.size() + 1; ++j) {
			int adj = 0;
			if (i == seq1.size() || j == seq2.size())
				adj = 1;

			lr = psa_mats.back()->at(i - 1)[j - 1].val + value_ref[pos][j - 1];			//diagonal (alignment)
			ud = lr;

			int gap_lr = 0, gap_ud = 0, pes_gap;

			for (int k = i - 1; k >= 0; --k) {
				double cand;
				if (psa_mats.back()->at(k)[j].gap >= ext_gap_ref.size() || i - k > ext_gap_ref[0].size()) {
					cand = psa_mats.back()->at(k)[j].val - ext_pen * (i - k);
					if (cand >= ud) {
						ud = cand;
						gap_ud = 1;
					}
				}
				else {
					cand = psa_mats.back()->at(k)[j].val + ext_gap_ref[psa_mats.back()->at(k)[j].gap][i - k];
					if (cand >= ud) {
						ud = cand;
						gap_ud = 1;
					}
				}

				cand = psa_mats.back()->at(k)[j - 1].val + value_ref[subsref[seq1[k]]][j - 1] + ext_gap_ref[1][i - k - 1];

				if (k == i - 1 || cand >= pes_val) {
					pes_val = cand;
					if (i - k - 1 > 0)
						pes_gap = 1;
					else
						pes_gap = 0;
				}
			}

			if (i < seq1.size())
				adj = 0;

			for (int k = j - 1; k >= 0; --k) {
				double cand;
				if (psa_mats.back()->at(i)[k].gap >= ext_gap_ref.size() || j - k > ext_gap_ref[0].size()) {
					cand = psa_mats.back()->at(i)[k].val - ext_pen * (j - k);
					if (cand >= lr) {
						lr = cand;
						gap_lr = 1;
					}
				}
				else {
					cand = psa_mats.back()->at(i)[k].val + ext_gap_ref[psa_mats.back()->at(i)[k].gap][j - k];
					if (cand >= lr) {
						lr = cand;
						gap_lr = psa_mats.back()->at(i)[k].gap + j - k;
					}
				}

				cand = psa_mats.back()->at(i)[k].pes_val - ext_pen * (j - k);

				if (cand >= pes_val) {
					pes_val = cand;
					pes_gap = 1;
				}
			}

			if (ud > lr)
				psa_mats.back()->at(i).emplace_back(ud, gap_ud, pes_val, pes_gap);
			else if (lr > ud || gap_lr >= gap_ud)
				psa_mats.back()->at(i).emplace_back(lr, gap_lr, pes_val, pes_gap);
			else
				psa_mats.back()->at(i).emplace_back(ud, gap_ud, pes_val, pes_gap);

		}
	}

	/*
	cout << "   -  ";
	for (int i = 0; i < seq2.size(); ++i)
	cout << (char)seq2[i] << "   ";
	cout << endl;
	for (int i = 0; i < psa_mats.back()->size(); ++i) {
	if (i > 0)
	cout << (char)seq1[i - 1] << " ";
	else
	cout << "   ";
	for (int j = 0; j < psa_mats.back()->at(i).size(); ++j) {
	cout << psa_mats.back()->at(i)[j].val << " ";
	}
	cout << endl;
	}
	cout << endl;
	for (int i = 0; i < psa_mats.back()->size(); ++i) {
	for (int j = 0; j < psa_mats.back()->at(i).size(); ++j) {
	cout << psa_mats.back()->at(i)[j].gap << " ";
	}
	cout << endl;
	}

	cout << "--------------------------------------\n";

	cout << "   -  ";
	for (int i = 0; i < seq2.size(); ++i)
	cout << (char)seq2[i] << "   ";
	cout << endl;
	for (int i = 0; i < psa_mats.back()->size(); ++i) {
	if (i > 0)
	cout << (char)seq1[i - 1] << " ";
	else
	cout << "   ";
	for (int j = 0; j < psa_mats.back()->at(i).size(); ++j) {
	cout << psa_mats.back()->at(i)[j].pes_val << " ";
	}
	cout << endl;
	}
	cout << endl;
	for (int i = 0; i < psa_mats.back()->size(); ++i) {
	for (int j = 0; j < psa_mats.back()->at(i).size(); ++j) {
	cout << psa_mats.back()->at(i)[j].pes_gap << " ";
	}
	cout << endl;
	}

	cout << "   -  ";
	for (int i = 0; i < seq2.size(); ++i)
	cout << (char)seq2[i] << "   ";
	cout << endl;
	for (int i = 0; i < psa_mats.back()->size(); ++i) {
	if (i > 0)
	cout << (char)seq1[i - 1] << " ";
	else
	cout << "   ";
	for (int j = 0; j < psa_mats.back()->at(i).size(); ++j) {
	cout << (psa_mats.back()->at(i)[j].pes_val < psa_mats.back()->at(i)[j].val) << " ";
	}
	cout << endl;
	}
	cout << endl;*/

	cout << "optimal psa value of " << psa_mats.size() << ": " << psa_mats.back()->back().back().val << endl;

	opt_psa += psa_mats.back()->back().back().val;

}
