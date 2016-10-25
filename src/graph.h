#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

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
