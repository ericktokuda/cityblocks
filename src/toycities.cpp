#include "toycities.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <set>
#include <cstdio>
#include <stdlib.h>
#include <igraph.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>

//export LD_PRELOAD=$HOME/.local/igraph-0.7.1/lib/libigraph.so
//clang++ src/toycities.cpp -I${HOME}/.local/igraph-0.7.1/include/igraph/ -L${HOME}/.local/igraph-0.7.1/lib/ -ligraph -o run


float mymean(vector<float> v) {// Not checking errors
	if (v.size() == 0)
	  throw "Mean of empty vector";
	float sum = 0;
	for(int i = 0; i < v.size(); i++)
		sum += v[i];
	return sum / v.size();
}

float mystd(vector<float> v, float avg) {// Not checking errors
	if (v.size() == 0)
	  throw "Std of empty vector";
	else if (v.size() == 1)
		return 0.0;

	float acc = 0;

	for(unsigned int i = 0; i < v.size(); i++)
		acc += pow(v[i] - avg, 2);
	return sqrt(acc / (v.size() - 1));
}

float mydiventropy(vector<float> v) {
	float acc = 0.0, sum_ = 0.0;
	for(unsigned int i=0; i < v.size(); i++)
		sum_ += v[i];

	for(unsigned int i=0; i < v.size(); i++) {
		float a = v[i] / sum_;
		acc -= a * log(a);
	}
	return acc;
}

float myevenness(vector<float> v) {
	float diventr = mydiventropy(v);
	return exp(diventr) / v.size();
}
// A utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
int compute_min_distance(vector<int> dist, vector<bool> sptSet) {
	// Initialize min value
	int min = INT_MAX, min_index;

	for (int v = 0; v < dist.size(); v++)
		if (sptSet[v] == false && dist[v] <= min)
			min = dist[v], min_index = v;

	return min_index;
}

// A utility function to print the constructed distance array
void printSolution(vector<int> dist) {
	printf("Vertex \t\t Distance from Source\n");
	for (int i = 0; i < dist.size(); i++)
		printf("%d \t\t %d\n", i, dist[i]);
}

// Function that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
float dijkstra(vector<vector<int>> graph, int src) {
	int n = graph.size();
	vector<int> dist(n, INT_MAX);
	vector<bool> sptSet(n, false);

	dist[src] = 0;

	for (int count = 0; count < n; count++) {
		int u = compute_min_distance(dist, sptSet);

		sptSet[u] = true;

		for (int v = 0; v < n; v++)
			if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX
				&& dist[u] + graph[u][v] < dist[v])
				dist[v] = dist[u] + graph[u][v];
	}

	//printSolution(dist);
	int acc = 0, nfinite = 0;
	for (int i = 0; i < dist.size(); i++) {
		if (dist[i] > 0 && dist[i] < INT_MAX) {
			acc += dist[i];
			nfinite ++;
		}
	}
	return (float) acc / nfinite;
}

float compute_average_path_length(vector<vector<int>> graph) {
	int n = graph.size();
	float acc = 0;

	for (int i = 0; i < n; i++)
		acc += dijkstra(graph, i);
	return acc / n;
}

//##########################################################

template <class T>
vector<T> concat(vector<T> v1, vector<T> v2) {
	v1.insert(v1.end(), v2.begin(), v2.end());
	return v1;
}
// Debugging function to print Fblock inner structures
void print_fblock(Fblock fblock) {
	printf("Fblockid: %d\nnodes:", fblock.id);
	for (int i = 0; i < fblock.nodes.size(); i++)
		printf("%d,", fblock.nodes[i]);
	printf(" neighbours:");
	for (int i = 0; i < fblock.neigh.size(); i++)
		printf("%d,", fblock.neigh[i]);
	printf(" edges:");
	for (int i = 0; i < fblock.edges.size(); i++)
		printf("%d,", fblock.edges[i]);
	printf("\n");
}

// Debugging function to print Block inner structures
void print_block(Block block) {
	printf("Blockid: %d\nfblocks:", block.id);
	for (int i = 0; i < block.fblocks.size(); i++)
		printf("%d,", block.fblocks[i]);
	printf("\n");
}

// Get the 4-conn. neighbours of node i, considering a regular grid
vector<int> get_4connected_neighbours(int i, int nrows, int ncols) {
	vector<int> neighbours;
	if (i >= ncols)
		neighbours.push_back(i - ncols);
	if ((i % ncols) != (ncols - 1))
		neighbours.push_back(i + 1);
	if (i < (nrows-1)*ncols)
		neighbours.push_back(i + ncols);
	if ((i % ncols) != 0)
		neighbours.push_back(i - 1);
	return neighbours;
}


// Get the nodes from a rectangular grid, given @nodesrows and @nodescols
// Returns a vector of Node structures
vector<Node> get_grid_nodes(int nodesrows, int nodescols) {
	vector<Node> nodes;
	for (int i = 0; i < nodesrows*nodescols; i++) {
		Node aux = {i , i / nodescols, i % nodescols};
		nodes.push_back(aux);
	}
	return nodes;
}


// Assuming the order defined in get_fundamental_blocks
vector<int> get_nodes_from_fblock(int fblockid,
								  int fblockrows,
								  int fblockcols) {
	if (fblockrows < 1 || fblockcols < 1)
		return vector<int> ();

	int topleft = fblockid + (fblockid / fblockcols);
	return vector<int> {topleft, topleft+1,
		topleft + 2 + fblockcols, topleft + 1 + fblockcols};
}


// Assuming the order defined in get_fundamental_blocks
vector<int> get_edges_from_fblock(int fblockid,
								  int fblockrows,
								  int fblockcols) {
	if (fblockrows < 1 || fblockcols < 1)
		return vector<int> ();

	//int nodesrows = fblockrows + 1, nodescols = fblockcols + 1;
	int top = fblockid;
	int nhorizedges = fblockcols * (fblockrows + 1);
	int left = (fblockid / fblockcols) + fblockid + nhorizedges;
	return vector<int> {top, left+1, top+fblockcols, left};
}



// Get the blocks defined by the grid
vector<Fblock> get_fundamental_blocks(int fblocksrows,
									  int fblockscols) {

	if (fblocksrows < 1 || fblockscols < 1)
		return vector<Fblock>();

	vector<Fblock> fblocks;

	for (int i = 0; i < fblocksrows; i++) {
		for (int j = 0; j < fblockscols; j++) {
			Fblock fblock;
			fblock.id = i*fblockscols+j;
			fblock.nodes = get_nodes_from_fblock(fblock.id,
												 fblocksrows, fblockscols);
			fblock.edges =  get_edges_from_fblock(fblock.id,
												  fblocksrows, fblockscols);
			int fblockpos = i * fblockscols + j;
			fblock.neigh = get_4connected_neighbours(fblockpos,
													 fblocksrows, fblockscols);

			fblocks.push_back(fblock);
		}
	}
	return fblocks;
}

//

// Assigns each Fblock structure to a single Block (1:1 relationship)
vector<Block> initialize_blocks(vector<Fblock> fblocks) {
	vector<Block> blocks;
	for (int i = 0; i < fblocks.size(); i++) {
		Fblock fblock = fblocks[i];
		Block block;
		block.id = i;
		block.fblocks.push_back(fblock.id);
		for (int j = 0; j < fblock.edges.size(); j++)
			block.edges.push_back(fblock.edges[j]);
		//block.neigh = fblocks[i].neigh;
		blocks.push_back(block);
	}
	return blocks;
}


// Assigns each Fblock structure to its Block container.
// Vector indices corresponds to the fblock ids.
vector<int> initialize_fblocks_ownership(vector<Block> blocks) {
	vector<int> ownership;
	map<int, int> ownershiphash;
	for (int i = 0; i < blocks.size(); i++) {
		int blockid = blocks[i].id;
		Block block = blocks[i];
		for (int j = 0; j < block.fblocks.size(); j++) {
			int fblockid = block.fblocks[j];
			//ownershiphash[fblockid] = block.fblocks[j];
			ownershiphash[fblockid] = blockid;
		}
	}

	// Assuming a continuous map of fblocks->block exists
	for (int i = 0; i < ownershiphash.size(); i++) {
		int blockid = ownershiphash[i];
		ownership.push_back(blockid);
	}
	return ownership;
}

// Get the edges of the rectangular grid (horiz. then vertical,
// from left-> right, top->bottom
// Order here is important because later we are gonna get the
// edge ids based on this order
vector<Edge> get_edges_from_regular_grid(int nodesrows, int nodescols) {
	int edgeid = 0;
	vector<Edge> edges;

	// First the horizontal edges
	for (int i = 0; i < nodesrows; i++) {
		for (int j = 0; j < nodescols-1; j++) {
			int leftnode = i*nodescols + j;
			Edge edg = {edgeid++, leftnode, leftnode + 1};
			edges.push_back(edg);
		}
	}

	// Then the vertical edges
	for (int i = 0; i < nodesrows-1; i++) {
		for (int j = 0; j < nodescols; j++) {
			int topid = i*nodescols + j;
			Edge edg = {edgeid++, topid, topid + nodescols};
			edges.push_back(edg);
		}
	}

	return edges;
}


vector<int> get_neighbour_blocks(Block block, vector<Fblock> fblocks,
								 vector<int> fblocksownership, int nrows, int ncols) {
	//vector<int> neighfblocks;
	vector<int> neighblocks;
	neighblocks.reserve(4*block.fblocks.size());

	for (int i = 0; i < block.fblocks.size(); i++) {
		vector<int> aux = get_4connected_neighbours(block.fblocks[i], nrows, ncols);

		//printf("neighblocks:");

		for (int j = 0; j < aux.size(); j++) {
			int neighid = fblocksownership[aux[j]];
			if (neighid != block.id)
				neighblocks.push_back(neighid);
			//printf("%d,", aux[j]);
		}
	}
	return neighblocks;
}


template <class T>
int get_idx_from_id(int id, T x) {
	for (int i = 0; i < x.size(); i++)
		if (x[i].id == id) return i;
	return -1;
}

template <class T>
vector<T> merge_vectors(vector<T> v1, vector<T> v2) {
	v1 = concat<T>(v1, v2);
	sort(v1.begin(), v1.end());
	v1.erase(unique(v1.begin(), v1.end()), v1.end());
	return v1;
}


vector<int> symmetric_diff_of_vectors(vector<int> v1, vector<int> v2) {
	set<int> s1(v1.begin(), v1.end());
	set<int> s2(v2.begin(), v2.end());
	set<int> s3 = s1;

	for (set<int, greater<int>>::iterator it=s2.begin(); it != s2.end(); ++it) {
		pair<set<int>::iterator, bool> ret;
		ret = s3.insert(*it);
		if (ret.second == 0) {// In case it was not possible to insert
			set<int, greater<int>>::iterator itaux = ret.first;
			s3.erase(itaux);
		}
	}
	vector<int> v(s3.begin(), s3.end());
	return v;
}
vector<Block> merge_blocks(int blockidx, int neighidx,
						   vector<Block> blocks, vector<Fblock> fblocks) {
	int keepidx = blockidx;
	int mergeidx = neighidx;

	if (neighidx < blockidx) { // Choose which block to keep
		keepidx = neighidx;
		mergeidx = blockidx;
	}

	// Merge fblocks
	blocks[keepidx].fblocks = concat(blocks[keepidx].fblocks,
									 blocks[mergeidx].fblocks);
	// Merge edges
	vector<int> e1 = blocks[keepidx].edges;
	vector<int> e2 = blocks[mergeidx].edges;
	blocks[keepidx].edges = symmetric_diff_of_vectors(e1, e2);

	blocks.erase(blocks.begin() + mergeidx);
	return blocks;
}


vector<vector<int>> get_adjmatrix_from_map(vector<Block> blocks, vector<Fblock> fblocks,
										   vector<Edge> edges, int fblockrows,
										   int fblockcols) {

	int nodesrows = fblockrows + 1, nodescols = fblockcols + 1;
	int nnodes = nodesrows * nodescols;
	set<int> alledgeids;

	vector<vector<int>> adj(nnodes, vector<int>(nnodes, 0));

	for (int i = 0; i < blocks.size(); i++) {
		Block block = blocks[i];
		for (int j = 0; j < block.edges.size(); j++) {
			Edge edge = edges[block.edges[j]];
			adj[edge[1]][edge[2]] = 1;
			adj[edge[2]][edge[1]] = 1;
			//printf("%d-%d ", edge[1], edge[2]);
			alledgeids.insert(edge[0]);
		}
		//totaledges += block.edges.size();
	}
	//printf("nedges:%ld\n", alledgeids.size());
	return adj;
}

vector<int> get_all_edges_flattened(vector<Block> blocks, vector<Fblock> fblocks,
									vector<Edge> edges, int fblockrows,
									int fblockcols) {
	int nodesrows = fblockrows + 1, nodescols = fblockcols + 1;
	int nedges = (nodescols-1)*nodesrows + (nodesrows-1)*nodescols;
	vector<int> edgesflat(nedges*2, 0);
	set<int> alledgeids;

	for (int i = 0; i < blocks.size(); i++) {
		Block block = blocks[i];
		for (int j = 0; j < block.edges.size(); j++) {
			Edge edge = edges[block.edges[j]];
			alledgeids.insert(edge[0]);
		}
	}

	int edgeind = 0;
	for (set<int>::iterator it = alledgeids.begin(); it != alledgeids.end(); ++it) {
		Edge edge = edges[*it];
		edgesflat[edgeind++] = edge[1];
		edgesflat[edgeind++] = edge[2];
	}
	//printf("nblocks:%d, theonedges:%d, edgeind:%d\n", blocks.size(), nedges, edgeind);
	return edgesflat;
}

//########################################################## TESTS
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
				cout << x[k] << ',';
			}
			cout << endl;
		}
	}
}

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
			printf("i:%d id:%d pos:(%d, %d)\n", k, nodes[k][0],
				   nodes[k][1], nodes[k][2]);
		}
	}
}

void test_get_nodes_from_fblock() {
	int fblockids[] = {1, 2, 3};
	int fblockrows[] = {0, 2, 4};
	int fblockcols[] = {0, 2, 4};

	printf("Started %s...\n", __func__);
	for (int j = 0; j < 3; j++) {
		//setbuf(stdout, NULL);
		printf("For a grid:(%d, %d) and fblockid:%d, nodes:",
			   fblockrows[j], fblockcols[j], fblockids[j]);
		vector<int> nodes = get_nodes_from_fblock(fblockids[j],
												  fblockrows[j],
												  fblockcols[j]);
		for (int k = 0; k < nodes.size(); k++) {
			printf("%d,", nodes[k]);
		}
	}
	printf("\nFinished %s\n", __func__);
}

void test_get_edges_from_fblock() {
	int fblockids[] = {1, 2, 3};
	int fblockrows[] = {0, 2, 3};
	int fblockcols[] = {0, 2, 4};

	printf("Started %s...\n", __func__);
	for (int j = 0; j < 3; j++) {
		//setbuf(stdout, NULL);
		printf("For a grid:(%d, %d) and fblockid:%d, edges:",
			   fblockrows[j], fblockcols[j], fblockids[j]);
		vector<int> edges = get_edges_from_fblock(fblockids[j],
												  fblockrows[j],
												  fblockcols[j]);
		for (int k = 0; k < edges.size(); k++) {
			printf("%d,", edges[k]);
		}
		printf("\n");
	}
	printf("\nFinished %s\n", __func__);
}

void test_get_fundamental_blocks() {
	int testnrows[] = {1, 2, 3};
	int testncols[] = {0, 2, 4};
	printf("%s...\n", __func__);
	for (int k = 0; k < 3; k++) {
		int nodesrows = testnrows[k];
		int nodescols = testncols[k];
		vector<Fblock> fblocks = get_fundamental_blocks(nodesrows,
														nodescols);
		printf("For a grid (%d, %d):\n",  nodesrows, nodescols);
		for (int i = 0; i < fblocks.size(); i++) {
			printf("%d id:%d nodes:", i, fblocks[i].id);
			Fblock fblock = fblocks[i];
			for (int j = 0; j < fblock.nodes.size(); j++)
				printf("%d,", fblock.nodes[j]);
			printf(" edges:");
			for (int j = 0; j < fblock.edges.size(); j++)
				printf("%d,", fblock.edges[j]);
			printf(" neigh:");
			for (int j = 0; j < fblock.neigh.size(); j++)
				printf("%d,", fblock.neigh[j]);
			printf("\n");
		}
	}
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
			Fblock fblock = fblocks[i];

			printf("id:%d, nodes:", fblock.id);
			for (int j = 0; j < fblock.nodes.size(); j++)
				printf("%d,", fblock.nodes[j]);

			printf(" Edges:");
			for (int j = 0; j < fblock.edges.size(); j++)
				printf("%d,", fblock.edges[j]);
			printf(" Neigh:");
			for (int j = 0; j < fblock.neigh.size(); j++)
				printf("%d,", fblock.neigh[j]);
			printf("\n");
		}

		vector<Block> blocks = initialize_blocks(fblocks);
		printf("Output (initialized blocks):\n");
		for (int i = 0; i < blocks.size(); i++) {
			Block block = blocks[i];

			printf("id:%d, fblock:", block.id);
			for (int j = 0; j < block.fblocks.size(); j++)
				printf("%d,", block.fblocks[j]);

			printf("\n");
		}
	}
}

void test_get_edges_from_regular_grid() {
	int testnrows[] = {1, 2, 3};
	int testncols[] = {0, 2, 4};

	printf("%s...\n", __func__);
	for (int k = 0; k < 3; k++) {
		int nodesrows = testnrows[k];
		int nodescols = testncols[k];
		vector<Edge> edges = get_edges_from_regular_grid(nodesrows,
														 nodescols);
		printf("For a grid (%d, %d)\n",  nodesrows, nodescols);

		for (int i = 0; i < edges.size(); i++) {
			//Edge edge = edges[i];
			//int edgeid = edges[i][0];
			//int edgeu = edges[i][1];
			//int edgev = edges[i][2];

			//printf("id:%d, nodes:", fblock.id);
			//printf("id:%d (%d, %d) ", edges[i][0],
			//edges[i][1],
			//edges[i][2]
			//);
			printf("id:%d (%d, %d) ", edges[i][0],
				   edges[i][1],
				   edges[i][2]
				  );
			printf(" ");
		}
		printf("\n");
	}
}

void test_get_neighbour_blocks() {
	// TODO
}

void test_igraph() {
	igraph_real_t avg_path;
	igraph_t graph;
	igraph_vector_t dimvector;
	igraph_vector_t edges;
	int i;

	igraph_vector_init(&dimvector, 2);
	VECTOR(dimvector)[0]=30;
	VECTOR(dimvector)[1]=30;
	igraph_lattice(&graph, &dimvector, 0, IGRAPH_UNDIRECTED, 0, 1);

	igraph_rng_seed(igraph_rng_default(), 42);
	igraph_vector_init(&edges, 20);
	for (i=0; i<igraph_vector_size(&edges); i++) {
		VECTOR(edges)[i] = rand() % (int)igraph_vcount(&graph);
	}

	igraph_average_path_length(&graph, &avg_path, IGRAPH_UNDIRECTED, 1);
	printf("Average path length (lattice):            %f\n", (double) avg_path);

	igraph_add_edges(&graph, &edges, 0);
	igraph_average_path_length(&graph, &avg_path, IGRAPH_UNDIRECTED, 1);
	printf("Average path length (randomized lattice): %f\n", (double) avg_path);

	igraph_vector_destroy(&dimvector);
	igraph_vector_destroy(&edges);
	igraph_destroy(&graph);
}

void test() {
	//int fblockrows = 50, fblockcols = 50;
	//int nodesrows = fblockrows + 1, nodescols = fblockcols + 1;
	test_get_4connected_neighbours();
	test_get_grid_nodes();
	test_get_edges_from_fblock();
	test_get_fundamental_blocks();
	test_initialize_blocks();
	test_get_edges_from_regular_grid();
	printf("END\n");
}
//########################################################## TESTS END

//##########################################################
//int main(int, char*[]) {
////test(); return 0;

////int fblockrows = 3, fblockcols = 4;
//int fblockrows = 20, fblockcols = 30;
//int nodesrows = fblockrows + 1, nodescols = fblockcols + 1;

////srand(0);
//srand(time(0));

//vector<Edge> edges = get_edges_from_regular_grid(nodesrows, nodescols);
//vector<Fblock> fblocks = get_fundamental_blocks(fblockrows, fblockcols);
//vector<Block> blocks = initialize_blocks(fblocks);
//vector<int> fblockownership = initialize_fblocks_ownership(blocks);
//vector<vector<int>> adj = get_adjmatrix_from_map(blocks, fblocks, edges,
//fblockrows, fblockcols);
////ofstream _stream;
//FILE *fh = fopen("/tmp/results.csv", "w");
//fprintf(fh, "nblocks,avgpathlength,blocksmean,blocksstd,blockscv,blocksdiventr,degreestd,degreesnonnullstd\n");

////fprintf(fh, "%ld,%g,%f,%f,%f,%f\n", blocks.size(), avg_path, blocksmean,
////blocksstd, blockscv, blocksdiventr, degreestd);
//setbuf(fh, NULL);
////_stream.open ("/tmp/result.txt");
//for (int i = 0; i < 5000; i++) {
//if (blocks.size() == 1) break;
//// sample a block
//int blockidx = rand() % blocks.size();
//Block block = blocks[blockidx];
////int blockid = block.id;

//// get its neighbour blocks
//vector<int> neighsrepeated = get_neighbour_blocks(block,
//fblocks,
//fblockownership,
//fblockrows,
//fblockcols);

//// sample a neighbour block, weighted by the num of neighbour fblocks
//int neighid = neighsrepeated[rand() % neighsrepeated.size()];
//int neighidx = get_idx_from_id<vector<Block>>(neighid, blocks);

//// Merge two blocks (update variables)
//blocks = merge_blocks(blockidx, neighidx, blocks, fblocks);

//// Update fblockownership
//for (int j = 0; j < blocks.size(); j++) {
//Block bl = blocks[j];
//for (int jj = 0; jj < bl.fblocks.size(); jj++) {
//fblockownership[bl.fblocks[jj]] = bl.id;
//}
//}

////vector<vector<int>> adj = get_adjmatrix_from_map(blocks, fblocks, edges,
////fblockrows, fblockcols);
////float l = compute_average_path_length(adj);

//igraph_vector_t ig_edges;
//vector<int> edgesflat = get_all_edges_flattened(blocks, fblocks, edges,
//fblockrows, fblockcols);
//igraph_real_t edgesigraph[edgesflat.size()];
//copy(edgesflat.begin(), edgesflat.end(), edgesigraph);
//igraph_vector_view(&ig_edges, edgesigraph, (igraph_real_t)edgesflat.size());
//igraph_t ig_graph;
//igraph_integer_t nn = nodesrows*nodescols;
////int ret = igraph_create(&ig_graph, &ig_edges, (igraph_integer_t) nn,
//igraph_create(&ig_graph, &ig_edges, (igraph_integer_t) nn,
//IGRAPH_UNDIRECTED);
//igraph_real_t avg_path;
//igraph_average_path_length(&ig_graph, &avg_path, IGRAPH_UNDIRECTED, 1);

////_stream << << avg_path << endl;
////printf("%f\n", mymean(blocks));
//vector<float> blockareas;
//float areassum = 0;
//blockareas.reserve(blocks.size());
//for (int k = 0; k < blocks.size(); k++) {
//blockareas.push_back(blocks[k].fblocks.size());
//areassum += blocks[k].fblocks.size();
//}

////for (int k = 0; k < blockareas.size(); k++)
////printf("%f ", blockareas[k]);
////printf("\n");

//float blocksmean = mymean(blockareas);
//float blocksstd = mystd(blockareas, blocksmean);
//float blockscv = blocksstd / blocksmean;
//float blocksdiventr =  mydiventropy(blockareas, areassum);
////float evenness =  myevenness(blockareas, areassum);
//igraph_vector_t ig_degrees;
//igraph_vector_init(&ig_degrees, 0);
//int ret = igraph_degree(&ig_graph, &ig_degrees, igraph_vss_all(), IGRAPH_ALL,
//false);


////vector<int> degrees(igraph_vector_size(&ig_degrees), 0);
//vector<float> degrees, degreesnonnull;
//degrees.reserve(igraph_vector_size(&ig_degrees));
//for (int k = 0; k < igraph_vector_size(&ig_degrees); k++) {
//degrees.push_back(static_cast<float>(VECTOR(ig_degrees)[k]));
//if (VECTOR(ig_degrees)[k] > 0)
//degreesnonnull.push_back(static_cast<float>(VECTOR(ig_degrees)[k]));
//}

//float degreestd = mystd(degrees, mymean(degrees));
//float degreesnonnullstd = mystd(degreesnonnull, mymean(degreesnonnull));
////printf("%f ", blockareas[k]);
////printf("\n");
////
////for (int k = 0; k < degrees.size(); k++)
////printf("%d ", degrees[k]);
////printf("\n");
////printf("\n");
////printf("\n");


////int igraph_degree(const igraph_t *graph, igraph_vector_t *res,
////const igraph_vs_t vids,
////igraph_neimode_t mode, igraph_bool_t loops);

////printf("%f %f %f %f %f\n", mean_, std_, cv_, diventr, degreestd);

//fprintf(fh, "%ld,%g,%f,%f,%f,%f,%f,%f\n", blocks.size(), avg_path, blocksmean,
//blocksstd, blockscv, blocksdiventr, degreestd, degreesnonnullstd);
////igraph_delete_edges(igraph_t *graph, igraph_es_t edges);
////printf("nblocks:%ld, avgpathlength:%f\n", blocks.size(), l);
////print_adj_matrix(adj);
//igraph_destroy(&ig_graph);
//}
//return 0;
//}

//template int get_idx_from_id<int>(int id, int x);
template int get_idx_from_id<vector<Block>>(int id, vector<Block> x);
