#pragma once

#include <vector>
//#include "catch2/catch.hpp"
#include <igraph.h>

#include <array>
//#include <stdlib.h>
//#include <igraph.h>

using namespace std;

#include <limits.h>
#include <stdio.h>


typedef array<int, 3> Node; // id, x, y
typedef array<int, 3> Edge; // id, x, y

typedef struct {
	int id;
	vector<int> nodes; // In a counter-clockwise order, from topleft
	vector<int> edges; // In a counter-clockwise order, from top
	vector<int> neigh;
} Fblock; // fundamental block

typedef struct {
	int id;
	vector<int> fblocks;
	vector<int> edges; // In a counter-clockwise order, from top
} Block; // Agglomerate of blocks

int compute_min_distance(vector<int> dist, vector<bool> sptSet);
void printSolution(vector<int> dist);
float dijkstra(vector<vector<int>> graph, int src);
float compute_average_path_length(vector<vector<int>> graph);
template <class T>
vector<T> concat(vector<T> v1, vector<T> v2);
void print_fblock(Fblock fblock);
void print_block(Block block);
vector<int> get_4connected_neighbours(int i, int nrows, int ncols);
//void test_get_4connected_neighbours();
vector<Node> get_grid_nodes(int nodesrows, int nodescols);
//void test_get_grid_nodes();
vector<int> get_nodes_from_fblock(int fblockid,
								  int fblockrows,
								  int fblockcols);

//void test_get_nodes_from_fblock();
vector<int> get_edges_from_fblock(int fblockid,
								  int fblockrows,
								  int fblockcols);


//void test_get_edges_from_fblock();
vector<Fblock> get_fundamental_blocks(int fblocksrows,
									  int fblockscols);
//void test_get_fundamental_blocks();
vector<Block> initialize_blocks(vector<Fblock> fblocks);
//void test_initialize_blocks();
vector<int> initialize_fblocks_ownership(vector<Block> blocks);
vector<Edge> get_edges_from_regular_grid(int nodesrows, int nodescols);
//void test_get_edges_from_regular_grid();
vector<int> get_neighbour_blocks(Block block, vector<Fblock> fblocks,
								 vector<int> fblocksownership,
								 int nrows, int ncols);
//void test_get_neighbour_blocks();
template <class T>
int get_idx_from_id(int id, T x);
template <class T>
vector<T> merge_vectors(vector<T> v1, vector<T> v2);
vector<int> symmetric_diff_of_vectors(vector<int> v1, vector<int> v2);
vector<Block> merge_blocks(int blockidx, int neighidx,
						   vector<Block> blocks, vector<Fblock> fblocks);
//void test_igraph();
vector<vector<int>> get_adjmatrix_from_map(vector<Block> blocks,
										   vector<Fblock> fblocks,
										   vector<Edge> edges,
										   int fblockrows,
										   int fblockcols);
vector<int> get_all_edges_flattened(vector<Block> blocks,
									vector<Fblock> fblocks,
									vector<Edge> edges, int fblockrows,
									int fblockcols);
