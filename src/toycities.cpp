#include <iostream>
#include <vector>
#include <igraph.h>

using namespace std;
//clang++ src/toycities.cpp -I/tmp/del/include/igraph/ -L/tmp/del/lib/ -ligraph && ./a.out

vector<int> get_4connected_neighbours(int i, int nrows, int ncols) {
	vector<int> neighbours;
	if (i >= ncols)
        neighbours.push_back(i - ncols);
	if (i % ncols != 0)
        neighbours.push_back(i - 1);
	if (i % ncols != ncols - 1)
        neighbours.push_back(i + 1);
	if (i < (nrows-1)*ncols)
        neighbours.push_back(i + ncols);
	return neighbours;
}

typedef struct {
	int id;
	int pos[2];
	//vector<int> neigh;
} Node;

typedef struct {
	int id;
	vector<int> nodes;
	vector<int> neigh;
} Fblock; // fundamental block

typedef struct {
	int id;
	vector<int> fblocks;
} Block; // Agglomerate of blocks

void test_get_4connected_neighbours() {
	int nrows = 3, ncols = 4;
	for (int i = 0; i < nrows*ncols; i++) {
		vector<int> x =  get_4connected_neighbours(i, nrows, ncols);

		for (int i = 0; i < x.size(); i++) {
			cout << x.at(i) << ',';
		}
		cout << endl;
	}
}

vector<Node> get_lattice_nodes(int nodesrows, int nodescols) {
	vector<Node> nodes;
	for (int i = 0; i < nodesrows*nodescols; i++) {
		Node aux;
		aux.id = i;
		aux.pos[0] = i / nodescols;
		aux.pos[1] = i % nodescols;
		nodes.push_back(aux);
	}
	return nodes;
}

vector<Fblock> get_fundamental_blocks(int nodesrows, int nodescols) {
	vector<Fblock> fblocks;
	for (int i = 0; i < nodesrows - 1; i++) {
		for (int j = 0; j < nodescols -1; j++) {
			Fblock fblock;
			fblock.id = i*nodescols+j;
			fblock.nodes.push_back(fblock.id);
			fblock.nodes.push_back(fblock.id + 1);
			fblock.nodes.push_back(fblock.id + nodescols);
			fblock.nodes.push_back(fblock.id + nodescols + 1);
			int fblockpos = i * (nodescols - 1) + j;
			fblock.neigh = get_4connected_neighbours(fblockpos, nodesrows-1, nodescols-1);
			fblocks.push_back(fblock);
		}
	}
	return fblocks;
}

void test_get_fundamental_blocks() {
	vector<Fblock> fblocks = get_fundamental_blocks(4, 5);

	for (int i = 0; i < fblocks.size(); i++) {
		cout << (fblocks.at(i)).id << ",nodes:";
		for (int j = 0; j < fblocks.at(j).nodes.size(); j++) {
			cout << (fblocks.at(i).nodes.at(j)) << ',';
		}
		cout << endl;
	}
}

vector<Block> initialize_blocks(vector<Fblock> fblocks) {
	vector<Block> blocks;
	for (int i = 0; i < fblocks.size(); i++) {
		Block aux;
		aux.id = i;
		aux.fblocks.push_back(fblocks.at(i).id);
		blocks.push_back(aux);
	}
	return blocks;
}

vector<Block> initialize_blocks(vector<Fblock> fblocks) {
	vector<Block> blocks;
	for (int i = 0; i < fblocks.size(); i++) {
		Block aux;
		aux.id = i;
		aux.fblocks.push_back(fblocks.at(i).id);
		blocks.push_back(aux);
	}
	return blocks;
}

int main(int, char*[]) {
	int nodesrows = 4, nodescols = 5;
	vector<Node> nodes = get_lattice_nodes(nodesrows, nodescols);
	vector<Fblock> fblocks = get_fundamental_blocks(nodesrows, nodescols);
	vector<Block> blocks = initialize_blocks(fblocks);


	//for (int i = 0; i < fblocks.size(); i++) {
		//cout << (fblocks.at(i)).id << ",nodes:";
		//for (int j = 0; j < fblocks.at(j).nodes.size(); j++) {
			//cout << (fblocks.at(i).nodes.at(j)) << ',';
		//}
		//cout << endl;
	//}
}
