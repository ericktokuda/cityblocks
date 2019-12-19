#include "catch2/catch.hpp"
#include "toycities.hpp"

using namespace std;

TEST_CASE("mymean", "[mymean]") {
	vector<float> v(1, 10);
	CHECK( mymean(v) == 10 );	//[10]
	v.push_back(0.0);
	CHECK( mymean(v) == 5 );	//[10, 0]
	v.push_back(-10.0);
	CHECK( mymean(v) == 0 );	//[10, 0, -10]
	vector<float> v2(1, 0.000000000001);
	CHECK( mymean(v2) == Approx(0.000000000001 ));	//[0.00000001]
}

TEST_CASE("mystd", "[mystd]") {
	vector<float> v(1, 10);
	CHECK( mystd(v, 10) == 0 );	//[10]
	v.push_back(0.0);
	CHECK( Approx(mystd(v, 5)) == 7.071 );	//[10, 0]
}

TEST_CASE("mydiventropy", "[mydiventropy]") {
	vector<float> v(1, 10);
	CHECK( mydiventropy(v) == Approx(0) );

	vector<float> v2(5, 10);
	CHECK( mydiventropy(v2) == Approx(log(5)) );

	vector<float> v3{1,1,8};
	CHECK( mydiventropy(v3) == Approx(-0.1*log(0.1)-0.1*log(0.1)-0.8*log(0.8)) );
}

TEST_CASE("myevenness", "[myevenness]") {
	vector<float> v(1, 10);
	CHECK( mydiventropy(v) == Approx(0) );

	vector<float> v2(5, 10);
	CHECK( mydiventropy(v2) == Approx(log(5)) );

	vector<float> v3{1,1,8};
	CHECK( mydiventropy(v3) == Approx(-0.1*log(0.1)-0.1*log(0.1)-0.8*log(0.8)) );
}

TEST_CASE("get_4connected_neighbours", "[get_4connected_neighbours]") {
	vector<int> v1{1,3}, v2{1,5,7,3};
	CHECK_THAT( get_4connected_neighbours(0, 3, 3), Catch::Equals(v1) );
	CHECK_THAT( get_4connected_neighbours(4, 3, 3), Catch::Equals(v2) );
}

TEST_CASE("get_grid_nodes", "[get_grid_nodes]") {
	vector<Node> v1{}, v2{{0,0,0}}, v3{{0,0,0}, {1,0,1}, {2,0,2}, {3,1,0}, {4,1,1}, {5,1,2}};
	CHECK_THAT( get_grid_nodes(0, 0), Catch::Equals(v1) );
	CHECK_THAT( get_grid_nodes(1, 1), Catch::Equals(v2) );
	CHECK_THAT( get_grid_nodes(2, 3), Catch::Equals(v3) );
}

//vector<int> get_nodes_from_fblock(int fblockid,
								  //int fblockrows,
								  //int fblockcols) {
								  
TEST_CASE("get_nodes_from_fblock", "[get_nodes_from_block]") {
	vector<int> v1{}, v2{1,2,5,4}, v3{3,4,7,6};
	CHECK_THAT( get_nodes_from_fblock(8, 0, 0), Catch::Equals(v1) );
	CHECK_THAT( get_nodes_from_fblock(1, 1, 2), Catch::Equals(v2) );
	CHECK_THAT( get_nodes_from_fblock(2, 2, 2), Catch::Equals(v3) );
}

TEST_CASE("get_edges_from_fblock", "[get_edges_from_fblock]") {
	vector<int> v1{}, v2{1,6,3,5}, v3{2,10,4,9};
	CHECK_THAT( get_edges_from_fblock(8, 0, 0), Catch::Equals(v1) );
	CHECK_THAT( get_edges_from_fblock(1, 1, 2), Catch::Equals(v2) );
	CHECK_THAT( get_edges_from_fblock(2, 2, 2), Catch::Equals(v3) );
}

//Not testing get_fundamental_block; simple function and very hard to test

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
