#include "mdd.hpp"
#include "build_mdd.hpp"

void Layer::popul_layer(int siz, int lvl, int mdd_count, bool last) {
	if (lvl == 0) {
		nodes.emplace_back(new vector<Node*>);
		nodes.back()->emplace_back(new Node(siz + 1, mdd_count));						//root node		
	}
	else if (lvl == 1) {
		nodes.emplace_back(new vector<Node*>);
		nodes.back()->reserve(siz + 1);
		for (size_t i = 0; i < siz; ++i)
			nodes.back()->emplace_back(new Node(i + 1, lvl, mdd_count));				//level one nodes
		nodes.back()->emplace_back(new Node(siz + 1, lvl, mdd_count));				//node with state = |S2| + 1					
	}
	else if (!last) {
		nodes.emplace_back(new vector<Node*>);
		nodes.back()->reserve(siz + 1);
		for (size_t i = 0; i < siz; ++i)
			nodes.back()->emplace_back(new Node(i + 1, siz, lvl, mdd_count));			//general nodes		
		nodes.back()->emplace_back(new Node(siz + 1, lvl, mdd_count));				//node with state = |S2| + 1	
		nodes.push_back(new vector<Node*>);
		nodes.back()->reserve(siz);
		for (size_t j = 0; j < siz; ++j)
			nodes.back()->emplace_back(new Node(j + 1, lvl, mdd_count));
	}
	else {
		nodes.emplace_back(new vector<Node*>);
		nodes.back()->reserve(siz + 1);
		for (size_t i = 0; i < siz; ++i)
			nodes.back()->emplace_back(new Node(i + 1, siz, lvl, mdd_count, last));	//general nodes		
		nodes.back()->emplace_back(new Node(siz + 1, lvl));				//node with state = |S2| + 1
		nodes.emplace_back(new vector<Node*>);
		nodes.back()->reserve(siz);
		for (size_t j = 0; j < siz; ++j)
			nodes.back()->emplace_back(new Node(j + 1, siz, lvl, mdd_count, last));
	}


}

void MDD::create_layer(int siz, int lvl, int mdd_count, bool last) {
	layers.push_back(new Layer(lvl));
	layers.back()->popul_layer(siz, lvl, mdd_count, last);
}


void MDD::create_arcs_r(vector<vector<Node*>*>& nodes1, vector<vector<Node*>*>& nodes2, int pos, size_t mdd_count) {

	size_t siz = strings[MDDs[mdd_count]->seq_to]->size();
	int siz2 = strings[MDDs[mdd_count]->seq_fr]->size();
	double bound = incumb - opt_psa + psa_mats[mdd_count]->back().back().val;

	if (ext_gap_ref[0][1] + psa_mats[mdd_count]->back()[siz2 - nodes2[0]->back()->lvl].pes_val + 0.0001 >= bound) {
		nodes2[0]->back()->max_to[0] = ext_gap_ref[0][1];
		nodes1[0]->at(0)->ogap_child = new Arc(ext_gap_ref[0][1], nodes2[0]->back(), nodes1[0]->at(0));			//opening gap arc (|S2| + 1 -> |S2| + 1)
		nodes2[0]->back()->ogap_parent = nodes1[0]->at(0)->ogap_child;
		nodes2[0]->back()->lng_arc = nodes1[0]->at(0)->ogap_child;
	}


	for (size_t j = 0; j < siz; ++j) {
		if (value_ref[pos][j] + ext_gap_ref[0][siz - j - 1] + psa_mats[mdd_count]->at(nodes2[0]->at(j)->state - 1)[siz2 - nodes2[0]->at(j)->lvl].val + 0.0001 >= bound) {
			nodes2[0]->at(j)->max_to[0] = value_ref[pos][j] + ext_gap_ref[0][siz - j - 1];
			nodes1[0]->at(0)->children.emplace_back(new Arc(value_ref[pos][j] + ext_gap_ref[0][siz - j - 1], nodes2[0]->at(j), nodes1[0]->at(0)));						//alignment arcs (implied gap cost + alignment reward)	
			nodes2[0]->at(j)->parents.emplace_back(nodes1[0]->at(0)->children.back());
			nodes2[0]->at(j)->lng_arc = nodes1[0]->at(0)->children.back();
		}
	}

}

void MDD::create_arcs(vector<vector<Node*>*>& nodes1, vector<vector<Node*>*>& nodes2, int lvl, int pos, size_t mdd_count) {
	
	size_t siz = strings[MDDs[mdd_count]->seq_to]->size();
	double bound = incumb - opt_psa + psa_mats[mdd_count]->back().back().val;
	int siz2 = strings[MDDs[mdd_count]->seq_fr]->size();

	if (nodes1[0]->back()->egap_parent == NULL && nodes1[0]->back()->ogap_parent == NULL) {
		nodes1[0]->back()->active = 0;
	}
	else{ 
		bool active = 0;
		if (nodes1[0]->back()->max_to[0] + ext_gap_ref[1][1] + psa_mats[mdd_count]->back()[siz2 - nodes2[0]->back()->lvl].pes_val + 0.0001 >= bound) {
			nodes2[0]->back()->max_to[0] = nodes1[0]->back()->max_to[0] + ext_gap_ref[1][1];
			nodes1[0]->back()->egap_child = new Arc(ext_gap_ref[1][1], nodes2[0]->back(), nodes1[0]->back());		//extension arc for node state = |S2| + 1
			nodes2[0]->back()->egap_parent = nodes1[0]->back()->egap_child;
			nodes2[0]->back()->lng_arc = nodes1[0]->back()->egap_child;
			active = 1;
		}
		for (size_t j = 0; j < siz; ++j) {														//alignment arcs for node state = |S2| + 1
			if (nodes1[0]->back()->max_to[0] + value_ref[pos][j] + ext_gap_ref[1][siz - j - 1] + psa_mats[mdd_count]->at(nodes2[0]->at(j)->state - 1)[siz2 - nodes2[0]->at(j)->lvl].val + 0.0001 >= bound) {
				nodes2[0]->at(j)->max_to[0] = nodes1[0]->back()->max_to[0] + value_ref[pos][j] + ext_gap_ref[1][siz - j - 1];
				nodes1[0]->back()->children.emplace_back(new Arc(value_ref[pos][j] + ext_gap_ref[1][siz - j - 1], nodes2[0]->at(j), nodes1[0]->back()));
				nodes2[0]->at(j)->parents.emplace_back(nodes1[0]->back()->children.back());
				nodes2[0]->at(j)->lng_arc = nodes1[0]->back()->children.back();
				active = 1;
			}
		}
		if (!active) {
			if (nodes1[0]->back()->ogap_parent != NULL) 
				nodes1[0]->back()->ogap_parent->reward[0] = -IloInfinity;
			if (nodes1[0]->back()->egap_parent != NULL) 
				nodes1[0]->back()->egap_parent->reward[0] = -IloInfinity;
			nodes1[0]->back()->max_fr[0] = -IloInfinity;
			nodes1[0]->back()->max_to[0] = -IloInfinity;
		}
	}
	int strt;
	if (lvl == 1)
		strt = 0;
	else
		strt = 1;
	for (int k = strt; k >= 0; --k) {
		for (int i = siz - 1; i >= 0; --i) {													//defining alignment arcs for nodes state > 2
			if (nodes1[k]->at(i)->ogap_parent == NULL && nodes1[k]->at(i)->egap_parent == NULL && nodes1[k]->at(i)->parents.empty()) {
				nodes1[k]->at(i)->active = 0;
			}
			else {
				bool active = 0;
				for (size_t j = 0; j < i; ++j) {
					if (nodes1[k]->at(i)->max_to[0] + value_ref[pos][j] + ext_gap_ref[k][i - j - 1] + psa_mats[mdd_count]->at(nodes2[0]->at(j)->state - 1)[siz2 - nodes2[0]->at(j)->lvl].val + 0.0001 >= bound) {
						nodes1[k]->at(i)->children.emplace_back(new Arc(value_ref[pos][j] + ext_gap_ref[k][i - j - 1], nodes2[0]->at(j), nodes1[k]->at(i)));
						if (nodes2[0]->at(j)->parents.empty() || nodes2[0]->at(j)->max_to[0] < nodes1[k]->at(i)->max_to[0] + value_ref[pos][j] + ext_gap_ref[k][i - j - 1]) {
							nodes2[0]->at(j)->max_to[0] = nodes1[k]->at(i)->max_to[0] + value_ref[pos][j] + ext_gap_ref[k][i - j - 1];
							nodes2[0]->at(j)->lng_arc = nodes1[k]->at(i)->children.back();
						}
						nodes2[0]->at(j)->parents.emplace_back(nodes1[k]->at(i)->children.back());
						active = 1;
					}
				}
				if (nodes1[k]->at(i)->max_to[0] + ext_gap_ref[k][1] + psa_mats[mdd_count]->at(nodes2[1]->at(i)->state - 1)[siz2 - nodes2[1]->at(i)->lvl].pes_val + 0.0001 >= bound) {
					if (k == 0) {
						nodes1[k]->at(i)->ogap_child = new Arc(ext_gap_ref[k][1], nodes2[1]->at(i), nodes1[k]->at(i));		//defining gap arc for nodes state > 2
						nodes2[1]->at(i)->ogap_parent = nodes1[k]->at(i)->ogap_child;
					}
					else {
						nodes1[k]->at(i)->egap_child = new Arc(ext_gap_ref[k][1], nodes2[1]->at(i), nodes1[k]->at(i));		//defining gap arc for nodes state > 2
						nodes2[1]->at(i)->egap_parent = nodes1[k]->at(i)->egap_child;
					}
					if (k == strt || nodes2[1]->at(i)->egap_parent == NULL || nodes2[1]->at(i)->max_to[0] < nodes1[k]->at(i)->max_to[0] + ext_gap_ref[k][1]) {
						nodes2[1]->at(i)->max_to[0] = nodes1[k]->at(i)->max_to[0] + ext_gap_ref[k][1];
						if (k == 0)
							nodes2[1]->at(i)->lng_arc = nodes1[k]->at(i)->ogap_child;
						else
							nodes2[1]->at(i)->lng_arc = nodes1[k]->at(i)->egap_child;
					}
					active = 1;
				}
				if (!active) {
					if (nodes1[k]->at(i)->ogap_parent != NULL) 
						nodes1[k]->at(i)->ogap_parent->reward[0] = -IloInfinity;
					if (nodes1[k]->at(i)->egap_parent != NULL)
						nodes1[k]->at(i)->egap_parent->reward[0] = -IloInfinity;
					for (int del = 0; del < nodes1[k]->at(i)->parents.size(); ++del)
						nodes1[k]->at(i)->parents[del]->reward[0] = -IloInfinity;
					nodes1[k]->at(i)->max_fr[0] = -IloInfinity;
					nodes1[k]->at(i)->max_to[0] = -IloInfinity;
				}
			}
		}
	}
}

void MDD::create_arcs_t(vector<vector<Node*>*>& nodes1, vector<vector<Node*>*>& nodes2, int lvl, int pos, size_t mdd_count) {

	size_t siz = strings[MDDs[mdd_count]->seq_to]->size();
	double bound = incumb - opt_psa + psa_mats[mdd_count]->back().back().val;
	int siz2 = strings[MDDs[mdd_count]->seq_fr]->size();

	if (nodes1[0]->back()->egap_parent == NULL) {
		nodes1[0]->back()->active = 0;
	}
	else {
		bool active = 0;
		if (nodes1[0]->back()->max_to[0] + ext_gap_ref[1][1 + siz] + 0.0001 >= bound) {
			nodes2[0]->back()->max_to[0] = nodes1[0]->back()->max_to[0] + ext_gap_ref[1][1 + siz];
			nodes1[0]->back()->egap_child = new Arc(ext_gap_ref[1][1 + siz], nodes2[0]->back(), nodes1[0]->back());		//extension arc for node state = |S2| + 1
			nodes2[0]->back()->egap_parent = nodes1[0]->back()->egap_child;
			nodes2[0]->back()->lng_arc = nodes1[0]->back()->egap_child;
			active = 1;
		}
		for (size_t j = 0; j < siz; ++j) {														//alignment arcs for node state = |S2| + 1
			if (nodes1[0]->back()->max_to[0] + value_ref[pos][j] + ext_gap_ref[1][siz - j - 1] + ext_gap_ref[0][j] + 0.0001 >= bound) {
				nodes2[0]->at(j)->max_to[0] = nodes1[0]->back()->max_to[0] + value_ref[pos][j] + ext_gap_ref[1][siz - j - 1] + ext_gap_ref[0][j];
				nodes1[0]->back()->children.emplace_back(new Arc(value_ref[pos][j] + ext_gap_ref[1][siz - j - 1] + ext_gap_ref[0][j], nodes2[0]->at(j), nodes1[0]->back()));
				nodes2[0]->at(j)->parents.emplace_back(nodes1[0]->back()->children.back());
				nodes2[0]->at(j)->lng_arc = nodes1[0]->back()->children.back();
				active = 1;
			}
		}
		if (!active) {
			nodes1[0]->back()->egap_parent->reward[0] = -IloInfinity;
			nodes1[0]->back()->max_fr[0] = -IloInfinity;
			nodes1[0]->back()->max_to[0] = -IloInfinity;
		}
	}

	for (int k = 1; k >= 0; --k) {							//defining alignment arcs for nodes state > 2
		for (int i = siz - 1; i >= 0; --i) {
			if (nodes1[k]->at(i)->ogap_parent == NULL && nodes1[k]->at(i)->egap_parent == NULL && nodes1[k]->at(i)->parents.empty()) {
				nodes1[k]->at(i)->active = 0;
			}
			else {
				bool active = 0;
				for (size_t j = 0; j < i; ++j) {
					if (nodes1[k]->at(i)->max_to[0] + value_ref[pos][j] + ext_gap_ref[k][i - j - 1] + ext_gap_ref[0][j] + 0.0001 >= bound) {
						nodes1[k]->at(i)->children.emplace_back(new Arc(value_ref[pos][j] + ext_gap_ref[k][i - j - 1] + ext_gap_ref[0][j], nodes2[0]->at(j), nodes1[k]->at(i)));
						if (nodes2[0]->at(j)->parents.empty() || nodes2[0]->at(j)->max_to[0] < nodes1[k]->at(i)->max_to[0] + value_ref[pos][j] + ext_gap_ref[k][i - j - 1] + ext_gap_ref[0][j]) {
							nodes2[0]->at(j)->max_to[0] = nodes1[k]->at(i)->max_to[0] + value_ref[pos][j] + ext_gap_ref[k][i - j - 1] + ext_gap_ref[0][j];
							nodes2[0]->at(j)->lng_arc = nodes1[k]->at(i)->children.back();
						}
						nodes2[0]->at(j)->parents.emplace_back(nodes1[k]->at(i)->children.back());
						active = 1;
					}
				}
				if (nodes1[k]->at(i)->max_to[0] + ext_gap_ref[k][i + 1] + 0.0001 >= bound) {
					if (k == 0) {
						nodes1[k]->at(i)->ogap_child = new Arc(ext_gap_ref[k][i + 1], nodes2[1]->at(i), nodes1[k]->at(i));		//defining gap arc for nodes state > 2
						nodes2[1]->at(i)->ogap_parent = nodes1[k]->at(i)->ogap_child;
					}
					else {
						nodes1[k]->at(i)->egap_child = new Arc(ext_gap_ref[k][i + 1], nodes2[1]->at(i), nodes1[k]->at(i));		//defining gap arc for nodes state > 2
						nodes2[1]->at(i)->egap_parent = nodes1[k]->at(i)->egap_child;
					}
					if (k == 1 || nodes2[1]->at(i)->egap_parent == NULL || nodes2[1]->at(i)->max_to[0] < nodes1[k]->at(i)->max_to[0] + ext_gap_ref[k][i + 1]) {
						nodes2[1]->at(i)->max_to[0] = nodes1[k]->at(i)->max_to[0] + ext_gap_ref[k][i + 1];
						if (k == 0)
							nodes2[1]->at(i)->lng_arc = nodes1[k]->at(i)->ogap_child;
						else
							nodes2[1]->at(i)->lng_arc = nodes1[k]->at(i)->egap_child;
					}
					active = 1;
				}
				if (!active) {
					if (nodes1[k]->at(i)->ogap_parent != NULL)
						nodes1[k]->at(i)->ogap_parent->reward[0] = -IloInfinity;
					if (nodes1[k]->at(i)->egap_parent != NULL)
						nodes1[k]->at(i)->egap_parent->reward[0] = -IloInfinity;
					for (int del = 0; del < nodes1[k]->at(i)->parents.size(); ++del)
						nodes1[k]->at(i)->parents[del]->reward[0] = -IloInfinity;
					nodes1[k]->at(i)->max_fr[0] = -IloInfinity;
					nodes1[k]->at(i)->max_to[0] = -IloInfinity;
				}
			}
		}
		/*if (nodes1[k]->empty()) {
			delete nodes1[k];
			nodes1.erase(nodes1.begin() + k);
		}
		else
			nodes1[k]->shrink_to_fit();*/
	}

	for (int k = nodes2.size() - 1; k >= 0; --k) {							//deleting last layer nodes with no parents
		for (int i = nodes2[k]->size() - 1; i >= 0; --i) {
			if (nodes2[k]->at(i)->ogap_parent == NULL && nodes2[k]->at(i)->egap_parent == NULL && nodes2[k]->at(i)->parents.empty())
				nodes2[k]->at(i)->active = 0;
		}
	}

	lng_path = layers.back()->nodes[0]->at(0);
	for (vector<vector<Node*>*>::iterator it = layers.back()->nodes.begin(); it != layers.back()->nodes.end(); ++it) {
		for (vector<Node*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2) 
			if ((*it2)->active && (*it2)->max_to[0] > lng_path->max_to[0])
				lng_path = (*it2);
	}
	cout << "%%%%%%%%%%%%%%%lng_tooo: " << lng_path->max_to[0] << endl;

}

void MDD::comp_lng_to(int rew_pos, bool restart) {

	for (size_t i = 1; i < layers.size(); ++i) {
		if (layers[i]->nodes.empty() || !layers[i]->active) {
			cout << "MDD is disconnected\n";
			return;
		}
		int aa = 0, bb = 0;
		for (vector<vector<Node*>*>::iterator it = layers[i]->nodes.begin(); it != layers[i]->nodes.end(); ++it) {
			//cout << "LAyer " << i << " node_layer" << aa << " ";
			for (vector<Node*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2) {
			//	cout << "Node " << bb << endl;
				if (restart || (*it2)->active)
					(*it2)->Comp_lng_to(rew_pos, restart);
				++bb;
			}
			++aa;
		}
	}

	//if (rew_pos == 0) {
		lng_path = layers.back()->nodes[0]->at(0);
		for (vector<vector<Node*>*>::iterator it = layers.back()->nodes.begin(); it != layers.back()->nodes.end(); ++it) {
			for (vector<Node*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
				if ((restart || (*it2)->active) && (*it2)->max_to[rew_pos] > lng_path->max_to[rew_pos])
					lng_path = (*it2);
		}
		cout << "lng_to: " << lng_path->max_to[rew_pos] << " ";
	//}

}


void MDD::comp_lng_fr(int rew_pos, bool restart) {

	for (int i = layers.size() - 2; i >= 0; --i) {
		if (layers[i]->nodes.empty() || !layers[i]->active) {
			cout << "MDD is disconnected\n";
			return;
		}
		int abc = 0;
		for (vector<vector<Node*>*>::iterator it = layers[i]->nodes.begin(); it != layers[i]->nodes.end(); ++it) {
			//cout << "gap level " << abc++ << endl;
			for (vector<Node*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
				if (restart || (*it2)->active)
					(*it2)->Comp_lng_fr(rew_pos, restart);
		}
	}

	cout << "lng_fr: " << layers[0]->nodes[0]->at(0)->max_fr[rew_pos] << " / ";

}

void MDD::extract_lngp() {

	Node* cur_nod = lng_path;

	for (int i = layers.size() - 1; i >= 1; --i) {
		vars[cur_nod->lng_arc->ID].val = 1;
		cur_nod = cur_nod->lng_arc->parent;
	}
}