#pragma once

#include <vector>
#include "load_inst.hpp"

class Arc;

class Node {
public:

	//size_t ID;
	int state;
	int lvl;
	int mdd_count;
	bool active;
	Arc* lng_arc;

	double max_to[2];				//solid arcs, opn gap, ext gaps ...
	double max_fr[2];

	Arc* gap_child;
	Arc* gap_parent;

	vector<Arc*> children;
	vector<Arc*> parents;

	void Comp_lng_to(int rew_pos, bool restart);
	void Comp_lng_fr(int rew_pos, bool restart);

	Node(int _state, int siz, int _lvl, int _mdd_count) {			//node constructor for general nodes
													
		state = _state;
		lvl = _lvl;
		active = 1;
		mdd_count = _mdd_count;

		gap_child = NULL;
		gap_parent = NULL;

		children.reserve(_state - 1);
		parents.reserve(siz - _state + 1);

	}

	Node(int _state, int siz, int _lvl, int _mdd_count, bool last) {			//node constructor for last layer nodes
			
		state = _state;
		lvl = _lvl;
		active = 1;
		mdd_count = _mdd_count;

		max_fr[0] = 0;
		max_fr[1] = 0;

		gap_child = NULL;
		gap_parent = NULL;

		parents.reserve(siz - _state + 1);
	}

	Node(int _state, int _mdd_count) {											//node constructor for root node
		state = _state;
		lvl = 0;
		active = 1;
		mdd_count = _mdd_count;

		max_to[0] = 0;
		max_to[1] = 0;

		gap_child = NULL;
		gap_parent = NULL;

		children.reserve(_state - 1);

	}

	Node(int _state, int _lvl, int _mdd_count) {								//node constructor for level one nodes and gap nodes (such nodes have one parent or one gap_parent)
		state = _state;
		lvl = _lvl;
		active = 1;
		mdd_count = _mdd_count;

		gap_child = NULL;
		gap_parent = NULL;

		children.reserve(abs(_state) - 1);

	}

	/*~Node() {
		for (vector<Arc*>::iterator it = children.begin(); it != children.end(); ++it)
			delete *it;
	}*/
};


class Arc {
public:
	double reward[2];
	Node* child;
	Node* parent;
	int ID;

	bool active;
	bool del_psa;

	Arc(double _rew, Node* chi, Node* par) { reward[0] = _rew; child = chi; parent = par; active = 1; del_psa = 0; ++num_arcs; }
};