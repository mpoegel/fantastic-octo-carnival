#include <cstdlib>
#include <iostream>

#include "graph.h"

using namespace std;


int main(int argc, char* argv[])
{

  if (argc < 2) {
    fprintf(stderr, "Usage: %s <input_graph>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char* graph_file = argv[1];

  int* srcs;
  int* dsts;
  int num_verts;
  int num_edges;
  int* out_array;
  int* out_degree_list;

  read_edge(graph_file, num_verts, num_edges, srcs, dsts);
  create_csr(num_verts, num_edges, srcs, dsts, out_array, out_degree_list);
  graph g = {num_verts, num_edges, out_array, out_degree_list};
  
  delete [] srcs;
  delete [] dsts;
  delete [] out_array;
  delete [] out_degree_list;

  return EXIT_SUCCESS;
}
