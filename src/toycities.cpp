#include <iostream>
#include <vector>
#include <igraph.h>

using namespace std;
//clang++ src/toycities.cpp -I/tmp/del/include/igraph/ -L/tmp/del/lib/ -ligraph && ./a.out

typedef struct {
	int id;
	int pos[2];
	//vector<int> neigh;
} Node;

typedef struct {
	int id;
	int uv[2];
	//vector<int> neigh;
} Edge;

typedef struct {
	int id;
	vector<int> nodes;
	vector<int> neigh;
} Fblock; // fundamental block

typedef struct {
	int id;
	vector<int> fblocks;
} Block; // Agglomerate of blocks

// Get the 4-connected neighbours of node i, considering a regular lattice
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

//
void test_get_4connected_neighbours() {
	int testnrows[] = {1, 2, 3};
	int testncols[] = {0, 2, 4};

	printf("%s...\n", __func__);
	for (int j = 0; j < 3; j++) {
		int nrows = testnrows[j];
		int ncols = testncols[j];
		printf("For a grid (%d, %d):\n",  nrows, ncols);
		for (int i = 0; i < nrows*ncols; i++) {
			cout << "Node " << i << ':';
			vector<int> x =  get_4connected_neighbours(i, nrows, ncols);
			for (int k = 0; k < x.size(); k++) {
				cout << x.at(k) << ',';
			}
			cout << endl;
		}
	}
}

// Get the nodes from a rectangular grid, given @nodesrows and @nodescols
// Returns a vector of Node structures
vector<Node> get_grid_nodes(int nodesrows, int nodescols) {
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

//
void test_get_grid_nodes() {
	int testnrows[] = {1, 2, 3};
	int testncols[] = {0, 2, 4};

	printf("%s...\n", __func__);
	for (int j = 0; j < 3; j++) {
		int nrows = testnrows[j];
		int ncols = testncols[j];
		vector<Node> nodes =  get_grid_nodes(nrows, ncols);
		printf("For a grid (%d, %d):\n",  nrows, ncols);
		for (int k = 0; k < nodes.size(); k++) {
			printf("i:%d id:%d pos:(%d, %d)\n", k, nodes.at(k).id, nodes.at(k).pos[0],
					nodes.at(k).pos[1]);
		}
	}
}

// Get the blocks defined by the grid. Returns a vector of Fblock structures
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

//
void test_get_fundamental_blocks() {
	int testnrows[] = {1, 2, 3};
	int testncols[] = {0, 2, 4};

	printf("%s...\n", __func__);
	for (int k = 0; k < 3; k++) {
		int nodesrows = testnrows[k];
		int nodescols = testncols[k];
		vector<Fblock> fblocks = get_fundamental_blocks(nodesrows, nodescols);
		printf("For a grid (%d, %d):\n",  nodesrows, nodescols);

		for (int i = 0; i < fblocks.size(); i++) {
			printf("%d id:%d nodes:", i, fblocks.at(i).id);
			Fblock fblock = fblocks.at(i);
			for (int j = 0; j < fblock.nodes.size(); j++)
				printf("%d,", fblock.nodes.at(j));
			printf(" neigh:");
			for (int j = 0; j < fblock.neigh.size(); j++)
				printf("%d,", fblock.neigh.at(j));
			printf("\n");
		}
	}
}

// Assigns each Block structure to a corresponding fblock for each of element
// in the vector of Fblock structures
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

// Get the edges of the rectangular grid
vector<Edge> get_edges_from_regular_grid(int nodesrows, int nodescols) {
	int edgeid = 0;
	vector<Edge> edges;

	for (int i = 0; i < nodesrows - 1; i++) {
		for (int j = 0; j < nodescols -1; j++) {
			Edge e1, e2;
			int nodeid = i*nodescols + j;
			e1.id = edgeid++;
			e1.uv[0] = nodeid;
			e1.uv[1] = nodeid + 1;
			edges.push_back(e1);
			e2.id = edgeid++;
			e2.uv[0] = nodeid;
			e2.uv[1] = nodeid + nodescols;
			edges.push_back(e2);
		}
	}

	for (int i = 0; i < nodesrows - 1; i++) {
			int nodeid = (i+1)*nodescols -1;
			Edge e1;
			e1.id = edgeid++;
			e1.uv[0] = nodeid;
			e1.uv[1] = nodeid + nodescols;
			edges.push_back(e1);
	}

	for (int j = 0; j < nodescols - 1; j++) {
			int nodeid = j*nodescols;
			Edge e1;
			e1.id = edgeid++;
			e1.uv[0] = nodeid;
			e1.uv[1] = nodeid + 1;
			edges.push_back(e1);
	}
	return edges;
}

int main(int, char*[]) {
	int nodesrows = 4, nodescols = 5;
	//test_get_4connected_neighbours();
	//test_get_grid_nodes();
	test_get_fundamental_blocks();
	//vector<Node> nodes = get_grid_nodes(nodesrows, nodescols);
	//vector<Fblock> fblocks = get_fundamental_blocks(nodesrows, nodescols);
	//vector<Block> blocks = initialize_blocks(fblocks);
	//vector<Edge> x = get_edges_from_regular_lattice(nodesrows, nodescols);

	//for (int i = 0; i < fblocks.size(); i++) {
		//cout << (fblocks.at(i)).id << ",nodes:";
		//for (int j = 0; j < fblocks.at(j).nodes.size(); j++) {
			//cout << (fblocks.at(i).nodes.at(j)) << ',';
		//}
		//cout << endl;
	//}
}
