#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <string>

#include "utils.h"

using namespace std;


#define out_degree(g, n) (g->out_degree_list[n+1] - g->out_degree_list[n])
#define out_vertices(g, n) &g->out_array[g->out_degree_list[n]]
#define out_weights(g, n) &g->out_weights[g->out_degree_list[n]]


struct graph {
  int num_verts;
  int num_edges;
  int* out_array;
  int* out_degree_list;
  float* latitudes;
  float* longitudes;
  double* out_weights;
};

void read_edge(char* filename, int &num_verts, int &num_edges, int* &srcs, int* &dsts);
void read_vert_latlong(char* filename, float* &latitudes, float* &longitudes);
void create_csr(int num_verts, int num_edges, int* srcs, int* dsts, float* latitudes, 
                float* longitudes, int* &out_array, double* &out_weights, int* &out_degree_list);

// shortest_paths(graph* g, int** dist, int** next);

// centrality algorithms
double* centrality_index(graph* g);
int* calc_betweenness(graph* g, double** dist, int** next);
double* calc_closeness(graph* g, double** dist);
double* calc_pagerank(graph* g, unsigned int num_iter, double damping_factor);

bool betweenness_comp(int* a, int* b);
void output_info(graph* g, int* betweenness);

int* match_by_population(graph* g, double* CI, int* populations, int num_cities);

#endif
