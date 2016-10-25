#include <cstdlib>
#include <iostream>

#include "graph.h"

using namespace std;


int* read_populations(char* filename, int* num_lines)
{
  ifstream infile;
  string line;
  infile.open(filename);
  
  getline(infile, line, ' ');
  *num_lines = atoi(line.c_str());
  getline(infile, line);
  int num_cities = atoi(line.c_str());

  int* populations = new int[num_cities];
  for (unsigned int i=0; i<num_cities; ++i) {
    populations[i] = -1;
  }
  
  int id;
  int pop;
  for (unsigned int i=0; i<*num_lines; ++i) {
    getline(infile, line, ' ');
    id = atoi(line.c_str());
    getline(infile, line);
    pop = atoi(line.c_str());
    populations[id] = pop;
  }

  infile.close();

  return populations;
}


float** read_latlong(char* filename)
{
  ifstream infile;
  string line;
  infile.open(filename);

  getline(infile, line);
  int num_cities = atoi(line.c_str());
  
  float** latlong = new float*[num_cities];
  for (unsigned int i=0; i<num_cities; ++i) {
    latlong[i] = new float[2];
    getline(infile, line, ' ');
    getline(infile, line, ' ');
    latlong[i][0] = atof(line.c_str());
    getline(infile, line);
    latlong[i][1] = atof(line.c_str());
  }

  infile.close();

  return latlong;
}

string* read_cities(char* filename)
{
  ifstream infile;
  string line;
  infile.open(filename);

  getline(infile, line);
  int num_cities = atoi(line.c_str());

  string* cities = new string[num_cities];
  for (unsigned int i=0; i<num_cities; ++i) {
    getline(infile, line, ' ');
    getline(infile, line);
    cities[i] = line;
  }

  infile.close();

  return cities;
}


int main(int argc, char* argv[])
{

  if (argc < 6) {
    fprintf(stderr, "Usage: %s <input_graph> <latlong_graph> <pop_file> <latlong_file>"
            "<cities_file>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char* graph_file = argv[1];
  char* latlong_graph_file = argv[2];
  char* pop_file = argv[3];
  char* latlong_file = argv[4];
  char* cities_file = argv[5];

  int num_cities;
  int* populations = read_populations(pop_file, &num_cities);
  float** latlong = read_latlong(latlong_file);
  string* cities = read_cities(cities_file);

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
  delete [] populations;
  for (unsigned int i=0; i<num_cities; ++i) {
    delete [] latlong[i];
  }
  delete [] latlong;
  delete [] cities;

  return EXIT_SUCCESS;
}
