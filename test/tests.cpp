#include "catch2/catch.hpp"

#include "toycities.hpp"
#include <iostream>
//#include <fstream>
//#include <algorithm>
//#include <map>
//#include <set>
//#include <cstdio>

//#include <stdlib.h>
#include <igraph.h>


//#include "toycities.hpp"

using namespace std;

TEST_CASE( "Factorials are computed", "[factorial]" ) {
    CHECK( 1 == 1 );
    CHECK( 3 == 1 );
    CHECK( 2 == 2 );
    CHECK( 2 == 2 );
    CHECK( 2 == 2 );
    CHECK( 2 == 2 );
}

TEST_CASE( "Bla", "[foo]" ) {
    CHECK( 1 == 1 );
    CHECK( 3 == 1 );
    CHECK( 2 == 2 );
}

void print_hello(const char* mystr) {
	printf("%s", mystr);
}

//void test_get_4connected_neighbours() {
	//int testnrows[] = {1, 2, 3};
	//int testncols[] = {0, 2, 4};

	//printf("%s...\n", __func__);
	//for (int j = 0; j < 3; j++) {
		//int nrows = testnrows[j];
		//int ncols = testncols[j];
		//printf("For a grid (%d, %d):\n",  nrows, ncols);
		//for (int i = 0; i < nrows*ncols; i++) {
			//cout << "Node " << i << ':';
			//vector<int> x =  get_4connected_neighbours(i, nrows, ncols);
			//for (int k = 0; k < x.size(); k++) {
				//cout << x[k] << ',';
			//}
			//cout << endl;
		//}
	//}
//}

//void test_get_grid_nodes() {
	//int testnrows[] = {1, 2, 3};
	//int testncols[] = {0, 2, 4};

	//printf("%s...\n", __func__);
	//for (int j = 0; j < 3; j++) {
		//int nrows = testnrows[j];
		//int ncols = testncols[j];
		//vector<Node> nodes =  get_grid_nodes(nrows, ncols);
		//printf("For a grid (%d, %d):\n",  nrows, ncols);
		//for (int k = 0; k < nodes.size(); k++) {
			//printf("i:%d id:%d pos:(%d, %d)\n", k, nodes[k][0],
				   //nodes[k][1], nodes[k][2]);
		//}
	//}
//}

//void test_get_nodes_from_fblock() {
	//int fblockids[] = {1, 2, 3};
	//int fblockrows[] = {0, 2, 4};
	//int fblockcols[] = {0, 2, 4};

	//printf("Started %s...\n", __func__);
	//for (int j = 0; j < 3; j++) {
		////setbuf(stdout, NULL);
		//printf("For a grid:(%d, %d) and fblockid:%d, nodes:",
			   //fblockrows[j], fblockcols[j], fblockids[j]);
		//vector<int> nodes = get_nodes_from_fblock(fblockids[j],
												  //fblockrows[j],
												  //fblockcols[j]);
		//for (int k = 0; k < nodes.size(); k++) {
			//printf("%d,", nodes[k]);
		//}
	//}
	//printf("\nFinished %s\n", __func__);
//}

//void test_get_edges_from_fblock() {
	//int fblockids[] = {1, 2, 3};
	//int fblockrows[] = {0, 2, 3};
	//int fblockcols[] = {0, 2, 4};

	//printf("Started %s...\n", __func__);
	//for (int j = 0; j < 3; j++) {
		////setbuf(stdout, NULL);
		//printf("For a grid:(%d, %d) and fblockid:%d, edges:",
			   //fblockrows[j], fblockcols[j], fblockids[j]);
		//vector<int> edges = get_edges_from_fblock(fblockids[j],
												  //fblockrows[j],
												  //fblockcols[j]);
		//for (int k = 0; k < edges.size(); k++) {
			//printf("%d,", edges[k]);
		//}
		//printf("\n");
	//}
	//printf("\nFinished %s\n", __func__);
//}

//void test_get_fundamental_blocks() {
	//int testnrows[] = {1, 2, 3};
	//int testncols[] = {0, 2, 4};
	//printf("%s...\n", __func__);
	//for (int k = 0; k < 3; k++) {
		//int nodesrows = testnrows[k];
		//int nodescols = testncols[k];
		//vector<Fblock> fblocks = get_fundamental_blocks(nodesrows,
														//nodescols);
		//printf("For a grid (%d, %d):\n",  nodesrows, nodescols);
		//for (int i = 0; i < fblocks.size(); i++) {
			//printf("%d id:%d nodes:", i, fblocks[i].id);
			//Fblock fblock = fblocks[i];
			//for (int j = 0; j < fblock.nodes.size(); j++)
				//printf("%d,", fblock.nodes[j]);
			//printf(" edges:");
			//for (int j = 0; j < fblock.edges.size(); j++)
				//printf("%d,", fblock.edges[j]);
			//printf(" neigh:");
			//for (int j = 0; j < fblock.neigh.size(); j++)
				//printf("%d,", fblock.neigh[j]);
			//printf("\n");
		//}
	//}
//}

//void test_initialize_blocks() {
	//int testnrows[] = {1, 2, 3};
	//int testncols[] = {0, 2, 4};
	//printf("%s...\n", __func__);
	//for (int k = 0; k < 3; k++) {
		//int nodesrows = testnrows[k];
		//int nodescols = testncols[k];
		//vector<Fblock> fblocks = get_fundamental_blocks(nodesrows, nodescols);
		//printf("Input (fundamental blocks):\n");

		//for (int i = 0; i < fblocks.size(); i++) {
			//Fblock fblock = fblocks[i];

			//printf("id:%d, nodes:", fblock.id);
			//for (int j = 0; j < fblock.nodes.size(); j++)
				//printf("%d,", fblock.nodes[j]);

			//printf(" Edges:");
			//for (int j = 0; j < fblock.edges.size(); j++)
				//printf("%d,", fblock.edges[j]);
			//printf(" Neigh:");
			//for (int j = 0; j < fblock.neigh.size(); j++)
				//printf("%d,", fblock.neigh[j]);
			//printf("\n");
		//}

		//vector<Block> blocks = initialize_blocks(fblocks);
		//printf("Output (initialized blocks):\n");
		//for (int i = 0; i < blocks.size(); i++) {
			//Block block = blocks[i];

			//printf("id:%d, fblock:", block.id);
			//for (int j = 0; j < block.fblocks.size(); j++)
				//printf("%d,", block.fblocks[j]);

			//printf("\n");
		//}
	//}
//}

//void test_get_edges_from_regular_grid() {
	//int testnrows[] = {1, 2, 3};
	//int testncols[] = {0, 2, 4};

	//printf("%s...\n", __func__);
	//for (int k = 0; k < 3; k++) {
		//int nodesrows = testnrows[k];
		//int nodescols = testncols[k];
		//vector<Edge> edges = get_edges_from_regular_grid(nodesrows,
														 //nodescols);
		//printf("For a grid (%d, %d)\n",  nodesrows, nodescols);

		//for (int i = 0; i < edges.size(); i++) {
			////Edge edge = edges[i];
			////int edgeid = edges[i][0];
			////int edgeu = edges[i][1];
			////int edgev = edges[i][2];

			////printf("id:%d, nodes:", fblock.id);
			////printf("id:%d (%d, %d) ", edges[i][0],
			////edges[i][1],
			////edges[i][2]
			////);
			//printf("id:%d (%d, %d) ", edges[i][0],
				   //edges[i][1],
				   //edges[i][2]
				  //);
			//printf(" ");
		//}
		//printf("\n");
	//}
//}

//void test_get_neighbour_blocks() {
	//// TODO
//}

//void test_igraph() {
	//igraph_real_t avg_path;
	//igraph_t graph;
	//igraph_vector_t dimvector;
	//igraph_vector_t edges;
	//int i;

	//igraph_vector_init(&dimvector, 2);
	//VECTOR(dimvector)[0]=30;
	//VECTOR(dimvector)[1]=30;
	//igraph_lattice(&graph, &dimvector, 0, IGRAPH_UNDIRECTED, 0, 1);

	//igraph_rng_seed(igraph_rng_default(), 42);
	//igraph_vector_init(&edges, 20);
	//for (i=0; i<igraph_vector_size(&edges); i++) {
		//VECTOR(edges)[i] = rand() % (int)igraph_vcount(&graph);
	//}

	//igraph_average_path_length(&graph, &avg_path, IGRAPH_UNDIRECTED, 1);
	//printf("Average path length (lattice):            %f\n", (double) avg_path);

	//igraph_add_edges(&graph, &edges, 0);
	//igraph_average_path_length(&graph, &avg_path, IGRAPH_UNDIRECTED, 1);
	//printf("Average path length (randomized lattice): %f\n", (double) avg_path);

	//igraph_vector_destroy(&dimvector);
	//igraph_vector_destroy(&edges);
	//igraph_destroy(&graph);
//}
