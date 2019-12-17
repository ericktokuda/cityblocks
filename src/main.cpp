//#include "toycities.hpp"

//#include <iostream>
//#include <fstream>
//#include <algorithm>
//#include <map>
//#include <set>
//#include <cstdio>

//#include <stdlib.h>
//#include <igraph.h>

//using namespace std;
////##########################################################
//int main(int, char*[]) {
	////int fblockrows = 3, fblockcols = 4;
	//int fblockrows = 50, fblockcols = 50;
	//int nodesrows = fblockrows + 1, nodescols = fblockcols + 1;
	////test_get_4connected_neighbours();
	////test_get_grid_nodes();
	////test_get_edges_from_fblock();
	////test_get_fundamental_blocks();
	////test_initialize_blocks();
	////test_get_edges_from_regular_grid();
	////printf("END\n");
	////return 0;

	////srand(0);
	//srand(time(0));

	//vector<Edge> edges = get_edges_from_regular_grid(nodesrows, nodescols);
	//vector<Fblock> fblocks = get_fundamental_blocks(fblockrows, fblockcols);
	//vector<Block> blocks = initialize_blocks(fblocks);
	//vector<int> fblockownership = initialize_fblocks_ownership(blocks);
	//vector<vector<int>> adj = get_adjmatrix_from_map(blocks, fblocks, edges,
													 //fblockrows, fblockcols);
	////ofstream _stream;
	//FILE *fh = fopen("/tmp/result.txt", "w");
	//fprintf(fh, "nblocks,avgpathlength\n");
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
		////int neighidx = get_idx_from_id<vector<Block>>(neighid, blocks);
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
		//igraph_destroy(&ig_graph);
		
		////_stream << << avg_path << endl;
		//fprintf(fh, "%ld,%g\n", blocks.size(), avg_path);
		////igraph_delete_edges(igraph_t *graph, igraph_es_t edges);
		////printf("nblocks:%ld, avgpathlength:%f\n", blocks.size(), l);
		////print_adj_matrix(adj);
	//}
	//return 0;
//}

