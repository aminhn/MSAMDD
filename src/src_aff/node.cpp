#include "load_inst.hpp"
#include "node.hpp"
#include <algorithm>
#include "build_mdd.hpp"


void Node::Comp_lng_to(int rew_pos, bool restart) {

	bool max_set = 0;

	max_to[rew_pos] = -IloInfinity;

	if (!parents.empty()) {
		size_t j = 0;
		while (j < parents.size()) {
			if (restart || parents[j]->active) {
				max_set = 1;
				max_to[rew_pos] = parents[j]->parent->max_to[rew_pos] + parents[j]->reward[rew_pos];
				if (rew_pos == 0)
					lng_arc = parents[j];
				break;
			}
			++j;
		}
		for (size_t i = j + 1; i < parents.size(); ++i) {
			if ((restart || parents[i]->active) && parents[i]->parent->max_to[rew_pos] + parents[i]->reward[rew_pos] > max_to[rew_pos]) {
				max_to[rew_pos] = parents[i]->parent->max_to[rew_pos] + parents[i]->reward[rew_pos];
				if (rew_pos == 0)
					lng_arc = parents[i];
			}
		}
	}


	if (ogap_parent && (restart || ogap_parent->active)) {
		if (!max_set || ogap_parent->parent->max_to[rew_pos] + ogap_parent->reward[rew_pos] > max_to[rew_pos]) {
			max_to[rew_pos] = ogap_parent->parent->max_to[rew_pos] + ogap_parent->reward[rew_pos];
			max_set = 1;
			if (rew_pos == 0)
				lng_arc = ogap_parent;
		}
	}

	if (egap_parent && (restart || egap_parent->active)) {
		if (!max_set || egap_parent->parent->max_to[rew_pos] + egap_parent->reward[rew_pos] > max_to[rew_pos]) {
			max_to[rew_pos] = egap_parent->parent->max_to[rew_pos] + egap_parent->reward[rew_pos];
			if (rew_pos == 0)
				lng_arc = egap_parent;
		}
	}

	//cout << max_to[rew_pos] << endl;

}


void Node::Comp_lng_fr(int rew_pos, bool restart) {

	bool max_set = 0;

	/*if (!gap_child) {
	cout << (gap_child == nullptr) << " " << children.size() << endl;
	cin.get();
	}*/

	max_fr[rew_pos] = -IloInfinity;

	if (!children.empty()) {
		size_t j = 0;
		while (j < children.size()) {
			if (restart || children[j]->active) {
				max_set = 1;
				max_fr[rew_pos] = children[j]->child->max_fr[rew_pos] + children[j]->reward[rew_pos];
				break;
			}
			++j;
		}
		for (size_t i = j + 1; i < children.size(); i++) {
			if ((restart || children[i]->active) && children[i]->child->max_fr[rew_pos] + children[i]->reward[rew_pos] > max_fr[rew_pos]) {
				max_fr[rew_pos] = children[i]->child->max_fr[rew_pos] + children[i]->reward[rew_pos];
			}
		}
	}

	if (ogap_child && (restart || ogap_child->active)) {
		if (!max_set || ogap_child->child->max_fr[rew_pos] + ogap_child->reward[rew_pos] > max_fr[rew_pos]) {
			max_fr[rew_pos] = ogap_child->child->max_fr[rew_pos] + ogap_child->reward[rew_pos];
			max_set = 1;
		}
	}

	if (egap_child && (restart || egap_child->active)) {
		if (!max_set || egap_child->child->max_fr[rew_pos] + egap_child->reward[rew_pos] > max_fr[rew_pos]) {
			max_fr[rew_pos] = egap_child->child->max_fr[rew_pos] + egap_child->reward[rew_pos];
		}
	}

	//int aaa = 2;
	//if (max_fr[0] > psa_mats[aaa]->at(state - 1)[9 - lvl].pes_val)
	//	cout << "state: " << state << " lvl " << lvl << " " << max_fr[0] << " " << psa_mats[aaa]->at(state - 1)[9 - lvl].val << " " << psa_mats[aaa]->at(state - 1)[9 - lvl].pes_val << " " << (max_fr[0] > psa_mats[aaa]->at(state - 1)[9 - lvl].val) <<  endl;

}



