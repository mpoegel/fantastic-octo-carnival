#include <algorithm>
#include <cstdlib>
#include <limits>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

#include "utils.h"

using namespace std;


#define out_degree(g, n) (g->out_degree_list[n+1] - g->out_degree_list[n])
#define out_vertices(g, n) &g->out_array[g->out_degree_list[n]]

struct graph {
  int num_verts;
  int num_edges;
  int* out_array;
  int* out_degree_list;
  float* latitudes;
  float* longitudes;
};

void read_edge(char* filename, int &num_verts, int &num_edges, int* &srcs, int* &dsts);
void read_vert_latlong(char* filename, float* &latitudes, float* &longitudes);
void create_csr(int num_verts, int num_edges, int* srcs, int* dsts, int* &out_array,
                int* &out_degree_list);

// centrality algorithms
double* centrality_index(graph* g);
int* calc_betweenness(graph* g);

bool betweenness_comp(int* a, int* b);
void output_info(graph* g, int* betweenness);

int* match_by_population(graph* g, double* CI, int* populations, int num_cities);
