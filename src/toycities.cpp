#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <igraph.h>

using namespace std;
//export LD_PRELOAD=$HOME/.local/igraph-0.7.1/lib/libigraph.so
//clang++ src/toycities.cpp -I${HOME}/.local/igraph-0.7.1/include/igraph/ -L${HOME}/.local/igraph-0.7.1/lib/ -ligraph -o run

#include <limits.h>
#include <stdio.h>

// Number of vertices in the graph
#define V 9

// A utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
int minDistance(int dist[], bool sptSet[])
{
    // Initialize min value
    int min = INT_MAX, min_index;

    for (int v = 0; v < V; v++)
        if (sptSet[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;

    return min_index;
}

// A utility function to print the constructed distance array
void printSolution(int dist[])
{
    printf("Vertex \t\t Distance from Source\n");
    for (int i = 0; i < V; i++)
        printf("%d \t\t %d\n", i, dist[i]);
}

// Function that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
void dijkstra(int graph[4][4], int src)
{
    int dist[V]; // The output array.  dist[i] will hold the shortest
    // distance from src to i

    bool sptSet[V]; // sptSet[i] will be true if vertex i is included in shortest
    // path tree or shortest distance from src to i is finalized

    // Initialize all distances as INFINITE and stpSet[] as false
    for (int i = 0; i < V; i++)
        dist[i] = INT_MAX, sptSet[i] = false;

    // Distance of source vertex from itself is always 0
    dist[src] = 0;

    // Find shortest path for all vertices
    for (int count = 0; count < V - 1; count++) {
        // Pick the minimum distance vertex from the set of vertices not
        // yet processed. u is always equal to src in the first iteration.
        int u = minDistance(dist, sptSet);

        // Mark the picked vertex as processed
        sptSet[u] = true;

        // Update dist value of the adjacent vertices of the picked vertex.
        for (int v = 0; v < V; v++)

            // Update dist[v] only if is not in sptSet, there is an edge from
            // u to v, and total weight of path from src to  v through u is
            // smaller than current value of dist[v]
            if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX
                && dist[u] + graph[u][v] < dist[v])
                dist[v] = dist[u] + graph[u][v];
    }

    // print the constructed distance array
    printSolution(dist);
}

//typedef struct {
	//int id;
	//int pos[2];
//} Node;

typedef int Node[3]; // id, x, y
typedef int Edge[3]; // id, u, v

//typedef struct {
	//int id;
	//int uv[2];
	////vector<int> neigh;
//} Edge;

typedef struct {
	int id;
	vector<int> nodes; // In a counter-clockwise order
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
		Node aux = {i , i / nodescols, i % nodescols};
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
			printf("i:%d id:%d pos:(%d, %d)\n", k, nodes.at(k)[0],
					nodes.at(k)[1], nodes.at(k)[2]);
		}
	}
}

// Considering a grid of sizes nrows, ncols
int get_top_left_node(int fblockid, int ncols) {
	return fblockid - fblockid / ncols;
}

// Get the blocks defined by the grid
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
			int topleft = get_top_left_node(fblock.id, fblockscols);
			fblock.nodes.push_back(topleft); // clockwise 
			fblock.nodes.push_back(topleft + 1);
			fblock.nodes.push_back(topleft + fblockscols + 1);
			fblock.nodes.push_back(topleft + fblockscols);
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
			//Edge e1, e2;
			int nodeid = i*nodescols + j;

			Edge e1 = {edgeid++, nodeid, nodeid + 1};
			edges.push_back(e1);

			Edge e2 = {edgeid++, nodeid, nodeid + 1};
			edges.push_back(e1);

			//e2.id = edgeid++;
			//e2.uv[0] = nodeid;
			//e2.uv[1] = nodeid + nodescols;
			//edges.push_back(e2);
		}
	}

	for (int i = 0; i < nodesrows - 1; i++) {
			int nodeid = (i+1)*nodescols -1;
			Edge e1 = {edgeid++, nodeid, nodeid + nodescols};
			edges.push_back(e1);
	}

	for (int j = 0; j < nodescols - 1; j++) {
			int nodeid = (nodesrows-1)*nodescols + j;
			Edge e1 = {edgeid++, nodeid, nodeid + 1};
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
		vector<Edge> edges = get_edges_from_regular_grid(nodesrows,
				nodescols);
		printf("For a grid (%d, %d)\n",  nodesrows, nodescols);

		for (int i = 0; i < edges.size(); i++) {
			Edge edge = edges.at(i);

			//printf("id:%d, nodes:", fblock.id);
			printf("id:%d (%d, %d) ", edges.at(i).id,
					edges.at(i)[1],
					edges.at(i)[2]
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

// Assuming a regular grid
vector<Edge> get_edges_from_fblock(Fblock fblock) {
	vector<int> nodes = fblock.nodes;
	for (int i = 0; i < 4; i++) {
		nodes.at(i)
	}
}

//vector<Edge> get_edges_from_block(Block block) {
	//vector<Fblock> fblocks;
	//for (int i = 0; i < fblocks.size(); i++) {
		//vector<int> v;
		//fblocks.at(i).nodes
	//}
//}

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

	   /* Let us create the example graph discussed above */
    //int graph[V][V] = { { 0, 4, 0, 0, 0, 0, 0, 8, 0 },
                        //{ 4, 0, 8, 0, 0, 0, 0, 11, 0 },
                        //{ 0, 8, 0, 7, 0, 4, 0, 0, 2 },
                        //{ 0, 0, 7, 0, 9, 14, 0, 0, 0 },
                        //{ 0, 0, 0, 9, 0, 10, 0, 0, 0 },
                        //{ 0, 0, 4, 14, 10, 0, 2, 0, 0 },
                        //{ 0, 0, 0, 0, 0, 2, 0, 1, 6 },
                        //{ 8, 11, 0, 0, 0, 0, 1, 0, 7 },
                        //{ 0, 0, 2, 0, 0, 0, 6, 7, 0 } };

	int graph[4][4] = {
		{ 0, 0, 0, 1 },
		{ 0, 0, 1, 1 },
		{ 0, 1, 0, 1 },
		{ 1, 1, 1, 0 }
	};

    dijkstra(graph, 2);

	//test_igraph();
	return 0;
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
