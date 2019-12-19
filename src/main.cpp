#include "toycities.hpp"
#include <igraph.h>

int main(int, char*[]) {
	//test(); return 0;
	
	//int fblockrows = 3, fblockcols = 4;
	int fblockrows = 20, fblockcols = 30;
	int nodesrows = fblockrows + 1, nodescols = fblockcols + 1;

	//srand(0);
	srand(time(0));

	vector<Edge> edges = get_edges_from_regular_grid(nodesrows, nodescols);
	vector<Fblock> fblocks = get_fundamental_blocks(fblockrows, fblockcols);
	vector<Block> blocks = initialize_blocks(fblocks);
	vector<int> fblockownership = initialize_fblocks_ownership(blocks);
	vector<vector<int>> adj = get_adjmatrix_from_map(blocks, fblocks, edges,
													 fblockrows, fblockcols);
	//ofstream _stream;
	FILE *fh = fopen("/tmp/results.csv", "w");
	fprintf(fh, "nblocks,avgpathlength,blocksmean,blocksstd,blockscv,blocksdiventr,blocksevenness,degreestd,degreesnonnullstd\n");

	setbuf(fh, NULL);
	//_stream.open ("/tmp/result.txt");
	for (int i = 0; i < 5000; i++) {
		if (blocks.size() == 1) break;
		// sample a block
		int blockidx = rand() % blocks.size();
		Block block = blocks[blockidx];
		//int blockid = block.id;

		// get its neighbour blocks
		vector<int> neighsrepeated = get_neighbour_blocks(block,
														  fblocks,
														  fblockownership,
														  fblockrows,
														  fblockcols);

		// sample a neighbour block, weighted by the num of neighbour fblocks
		int neighid = neighsrepeated[rand() % neighsrepeated.size()];
		int neighidx = get_idx_from_id<vector<Block>>(neighid, blocks);

		// Merge two blocks (update variables)
		blocks = merge_blocks(blockidx, neighidx, blocks, fblocks);

		// Update fblockownership
		for (int j = 0; j < blocks.size(); j++) {
			Block bl = blocks[j];
			for (int jj = 0; jj < bl.fblocks.size(); jj++) {
				fblockownership[bl.fblocks[jj]] = bl.id;
			}
		}

		igraph_vector_t ig_edges;
		vector<int> edgesflat = get_all_edges_flattened(blocks, fblocks, edges,
														fblockrows, fblockcols);
		igraph_real_t edgesigraph[edgesflat.size()];
		copy(edgesflat.begin(), edgesflat.end(), edgesigraph);
		igraph_vector_view(&ig_edges, edgesigraph, (igraph_real_t)edgesflat.size());
		igraph_t ig_graph;
		igraph_integer_t nn = nodesrows*nodescols;
		//int ret = igraph_create(&ig_graph, &ig_edges, (igraph_integer_t) nn,
		igraph_create(&ig_graph, &ig_edges, (igraph_integer_t) nn,
								IGRAPH_UNDIRECTED);
		igraph_real_t avg_path;
		igraph_average_path_length(&ig_graph, &avg_path, IGRAPH_UNDIRECTED, 1);
		
		vector<float> blockareas;
		float areassum = 0;
		blockareas.reserve(blocks.size());
		for (int k = 0; k < blocks.size(); k++) {
			blockareas.push_back(blocks[k].fblocks.size());
			areassum += blocks[k].fblocks.size();
		}

		float blocksmean = mymean(blockareas);
		float blocksstd = mystd(blockareas, blocksmean);
		float blockscv = blocksstd / blocksmean;
		float blocksdiventr =  mydiventropy(blockareas);
		float blocksevenness =  myevenness(blockareas);
		igraph_vector_t ig_degrees;
		igraph_vector_init(&ig_degrees, 0);
		igraph_degree(&ig_graph, &ig_degrees, igraph_vss_all(), IGRAPH_ALL, false);
		
		vector<float> degrees, degreesnonnull; // Convert ig_degrees to vector
		degrees.reserve(igraph_vector_size(&ig_degrees));
		for (int k = 0; k < igraph_vector_size(&ig_degrees); k++) {
			degrees.push_back(static_cast<float>(VECTOR(ig_degrees)[k]));
			if (VECTOR(ig_degrees)[k] > 0)
				degreesnonnull.push_back(static_cast<float>(VECTOR(ig_degrees)[k]));
		}

		float degreestd = mystd(degrees, mymean(degrees));
		float degreesnonnullstd = mystd(degreesnonnull, mymean(degreesnonnull));

		fprintf(fh, "%ld,%g,%f,%f,%f,%f,%f,%f,%f\n", blocks.size(), avg_path,
				blocksmean, blocksstd, blockscv, blocksdiventr, blocksevenness,
				degreestd, degreesnonnullstd);

		igraph_destroy(&ig_graph);
	}
	return 0;
}


