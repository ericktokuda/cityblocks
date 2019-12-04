#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <igraph.h>

using namespace std;
//clang++ src/toycities.cpp -I/tmp/del/include/igraph/ -L/tmp/del/lib/ -ligraph && ./a.out

typedef struct {
	int id;
	int pos[2];
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
	//vector<int> neigh;
} Block; // Agglomerate of blocks


void print_fblock(Fblock fblock) {
	printf("Fblockid: %d\nnodes:", fblock.id);
	for (int i = 0; i < fblock.nodes.size(); i++)
		printf("%d,", fblock.nodes.at(i));
	printf(" neighbours:");
	for (int i = 0; i < fblock.neigh.size(); i++)
		printf("%d,", fblock.neigh.at(i));
	printf("\n");
}

void print_block(Block block) {
	printf("Blockid: %d\nfblocks:", block.id);
	for (int i = 0; i < block.fblocks.size(); i++)
		printf("%d,", block.fblocks.at(i));
	printf("\n");
}
// Get the 4-connected neighbours of node i, considering a regular grid
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
vector<Fblock> get_fundamental_blocks(int fblocksrows, int fblockscols) {
	if (fblocksrows < 1 || fblockscols < 1) {
		vector<Fblock> empty;
		return empty;
	}

	vector<Fblock> fblocks;
	for (int i = 0; i < fblocksrows; i++) {
		for (int j = 0; j < fblockscols; j++) {
			Fblock fblock;
			fblock.id = i*fblockscols+j;
			fblock.nodes.push_back(fblock.id);
			fblock.nodes.push_back(fblock.id + 1);
			fblock.nodes.push_back(fblock.id + fblockscols);
			fblock.nodes.push_back(fblock.id + fblockscols + 1);
			int fblockpos = i * (fblockscols - 1) + j;
			fblock.neigh = get_4connected_neighbours(fblockpos, fblocksrows,
					fblockscols);
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
		//aux.neigh = fblocks.at(i).neigh;
		blocks.push_back(aux);
	}
	return blocks;
}

void test_initialize_blocks() {
	int testnrows[] = {1, 2, 3};
	int testncols[] = {0, 2, 4};

	printf("%s...\n", __func__);
	for (int k = 0; k < 3; k++) {
		int nodesrows = testnrows[k];
		int nodescols = testncols[k];
		vector<Fblock> fblocks = get_fundamental_blocks(nodesrows, nodescols);
		printf("Input (fundamental blocks):\n");

		for (int i = 0; i < fblocks.size(); i++) {
			Fblock fblock = fblocks.at(i);

			printf("id:%d, nodes:", fblock.id);
			for (int j = 0; j < fblock.nodes.size(); j++)
				printf("%d,", fblock.nodes.at(j));

			printf(" Neigh:");
			for (int j = 0; j < fblock.neigh.size(); j++)
				printf("%d,", fblock.neigh.at(j));
			printf("\n");
		}

		vector<Block> blocks = initialize_blocks(fblocks);
		printf("Output (derived blocks):\n");
		for (int i = 0; i < blocks.size(); i++) {
			Block block = blocks.at(i);

			printf("id:%d, nodes:", block.id);
			for (int j = 0; j < block.fblocks.size(); j++)
				printf("%d,", block.fblocks.at(j));

			printf("\n");
		}
	}
}

// Assigns each Block structure to a corresponding fblock for each of element
// in the vector of Fblock structures
vector<int> initialize_fblocks_ownership(vector<Fblock> fblocks) {
	vector<int> ownership;
	for (int i = 0; i < fblocks.size(); i++)
		ownership.push_back(i);
	return ownership;
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
			int nodeid = (nodesrows-1)*nodescols + j;
			Edge e1;
			e1.id = edgeid++;
			e1.uv[0] = nodeid;
			e1.uv[1] = nodeid + 1;
			edges.push_back(e1);
	}
	return edges;
}

void test_get_edges_from_regular_grid() {
	int testnrows[] = {1, 2, 3};
	int testncols[] = {0, 2, 4};

	printf("%s...\n", __func__);
	for (int k = 0; k < 3; k++) {
		int nodesrows = testnrows[k];
		int nodescols = testncols[k];
		vector<Edge> edges = get_edges_from_regular_grid(nodesrows, nodescols);
		printf("For a grid (%d, %d)\n",  nodesrows, nodescols);

		for (int i = 0; i < edges.size(); i++) {
			Edge edge = edges.at(i);

			//printf("id:%d, nodes:", fblock.id);
			printf("id:%d (%d, %d) ", edges.at(i).id,
					edges.at(i).uv[0],
					edges.at(i).uv[1]
					);
			printf(" ");
		}
		printf("\n");
	}
}

vector<int> get_neighbour_blocks(Block block, vector<Fblock> fblocks,
		vector<int> fblocksownership, int nrows, int ncols) {
	//vector<int> neighfblocks;
	vector<int> neighblocks;
	neighblocks.reserve(4*block.fblocks.size());

	for (int i = 0; i < block.fblocks.size(); i++) {
		vector<int> aux = get_4connected_neighbours(block.fblocks.at(i), nrows, ncols);

		//printf("neighblocks:");

		for (int j = 0; j < aux.size(); j++) {
			int neighid = fblocksownership.at(aux.at(j));
			if (neighid != block.id)
				neighblocks.push_back(neighid);
			//printf("%d,", aux.at(j));
		}
	}
	return neighblocks;
}

void test_get_neighbour_blocks() {
	// TODO
}

template <class Vectortype>
int get_idx_from_id(int id, Vectortype x) {
	for (int i = 0; i < x.size(); i++)
		if (x.at(i).id == id) return i;
	return -1;
}

vector<Block> merge_blocks(int blockidx, int neighidx, vector<Block> blocks) {
	int keepidx = blockidx;
	int mergeidx = neighidx;

	if (neighidx < blockidx) {
		keepidx = neighidx;
		mergeidx = blockidx;
	}
	
	//printf("\nblockidx:%d, neighidx:%d\n", blockidx, neighidx);
	//printf("\nkeepidx:%d, mergeidx:%d\n", keepidx, mergeidx);
	//print_block(blocks.at(keepidx));
	//printf("checkpoint1\n");
	//print_block(blocks.at(mergeidx));
	//printf("checkpoint2\n");
	//print_block(blocks.at(keepidx));
	// Merge  fblocks
	blocks.at(keepidx).fblocks.insert(blocks.at(keepidx).fblocks.end(),
			blocks.at(mergeidx).fblocks.begin(),
			blocks.at(mergeidx).fblocks.end());
	//print_block(blocks.at(keepidx));

	// Merge neighbours
	//vector<int> neigh = blocks.at(keepidx).neigh;
	//neigh.insert(blocks.at(keepidx).neigh.end(), // concatenate
			//blocks.at(mergeidx).neigh.begin(),
			//blocks.at(mergeidx).neigh.end());


	//sort(neigh.begin(), neigh.end()); // remove duplicates
	//neigh.erase(unique(neigh.begin(), neigh.end( ), neigh.end()));


	//blocks.at(keepidx).fblocks.insert(fblocks1.end(), fblocks2.begin(), fblocks2.end());
			

	//printf("checkpoint10\n");
	//printf("\nMiddle %ld", blocks.at(keepidx).fblocks.size());
	blocks.erase(blocks.begin() + mergeidx);
	//printf("\nAfter %ld", blocks.at(keepidx).fblocks.size());
	//printf("checkpoint12\n");
	return blocks;
}

int main(int, char*[]) {
	int blocksrows = 20, blockscols = 30;
	//test_get_4connected_neighbours();
	//test_get_grid_nodes();
	//test_get_fundamental_blocks();
	//test_initialize_blocks();
	//test_get_edges_from_regular_grid();
	
	//srand(0);
	srand(time(0));
	
	vector<Node> nodes = get_grid_nodes(blocksrows-1, blockscols-1);
	vector<Fblock> fblocks = get_fundamental_blocks(blocksrows, blockscols);
	vector<Block> blocks = initialize_blocks(fblocks);
	vector<int> fblockownership = initialize_fblocks_ownership(fblocks);

	for (int i = 0; i < 5000; i++) {
		if (blocks.size() == 1) break;
		// sample a block
		int blockidx = rand() % blocks.size();
		Block block = blocks.at(blockidx);
		int blockid = block.id;


		// get its neighbour blocks
		vector<int> neighsrepeated = get_neighbour_blocks(block, fblocks,
				fblockownership, blocksrows, blockscols);

		// sample a neighbour block, weighted by the num of neighbour fblocks
		int neighid = neighsrepeated[rand() % neighsrepeated.size()];
		int neighidx = get_idx_from_id<vector<Block>>(neighid, blocks);
		//printf("neighsrepeated.size:%ld neighid(x):%d,%d ", neighsrepeated.size(),
				//neighid, neighidx);

		// Merge two blocks (update variables)
		blocks = merge_blocks(blockidx, neighidx, blocks);
		
		printf("After blocksz:%ld, block.fblocks.sz:%ld,", blocks.size(),
				blocks.at(blockidx <= neighidx ? blockidx : neighidx).fblocks.size());
		printf("\n");

		// Update fblockownership
		for (int j = 0; j < blocks.size(); j++) {
			Block bl = blocks.at(j);
			for (int jj = 0; jj < bl.fblocks.size(); jj++) {
				fblockownership.at(bl.fblocks.at(jj)) = bl.id;
			}
		}
	}
	return 0;
}
