#pragma once

#include "load_inst.hpp"
#include <vector>
#include "node.hpp"

class Layer {
public:

	bool active;
	vector<vector<Node*>*> nodes;

	void popul_layer(int siz, int lvl, int mdd_count, bool last);

	Layer(int lvl) {
		nodes.reserve(lvl + 1);
		active = 1;
	}

/*	~Layer() {
		for (vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
			delete *it;

	}*/


};


class MDD {
public:

	vector<Layer*> layers;
	int seq_fr;
	int seq_to;

	Node* lng_path;

	void comp_lng_to(int rew_pos, bool restart);
	void comp_lng_fr(int rew_pos, bool restart);
	void extract_lngp();

	void create_layer(int siz, int lvl, int mdd_count, bool last);

	void create_arcs(vector<vector<Node*>*>& nodes1, vector<vector<Node*>*>& nodes2, int lvl, int pos, size_t mdd_count);
	void create_arcs_r(vector<vector<Node*>*>& nodes1, vector<vector<Node*>*>& nodes2, int pos, size_t mdd_count);
	void create_arcs_t(vector<vector<Node*>*>& nodes1, vector<vector<Node*>*>& nodes2, int lvl, int pos, size_t mdd_count);

	MDD(int siz, int i, int j) {
		seq_fr = i;
		seq_to = j;
		layers.reserve(siz + 1);
	}

	~MDD() {
		for (vector<Layer*>::iterator it = layers.begin(); it != layers.end(); it++) {
			(*it)->~Layer();
			delete *it;
		}
	}

};