#include <algorithm>
#include <cstdlib>
#include <iostream>

#include "graph.h"
#include "utils.h"

using namespace std;


int* read_populations(char* filename, int* num_cities)
{
  ifstream infile;
  string line;
  infile.open(filename);
  
  getline(infile, line, ' ');
  int num_lines = atoi(line.c_str());
  getline(infile, line);
  *num_cities = atoi(line.c_str());

  int* populations = new int[*num_cities];
  for (unsigned int i=0; i<*num_cities; ++i) {
    populations[i] = -1;
  }
  
  int id;
  int pop;
  for (unsigned int i=0; i<num_lines; ++i) {
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


void calculate_error(graph* g, int* y_hat, float** latlong, int num_cities)
{
  double* dists = new double[num_cities];
  double total_dist = 0.0;
  unsigned int missing = 0;
  for (unsigned int i=0; i<num_cities; ++i) {
    int v = y_hat[i];
    if (v > 0) {
      double d = measureLatLongDist(g->latitudes[v], g->longitudes[v], latlong[i][0],
                                    latlong[i][1]);
      dists[i] = d;
      total_dist += d;
    } else {
      ++missing;
    }
  }
  sort(dists, dists + num_cities);
  double min_dist = dists[0];
  double max_dist = dists[num_cities - 1];
  double avg_dist = total_dist / (double)num_cities;
  double med_dist;
  int mid = num_cities / 2;
  if (num_cities % 2 == 0) {
    med_dist = (dists[mid] + dists[mid - 1]) / 2.0;
  } else {
    med_dist = dists[mid];
  }
  
  printf("Analysis on %d/%d cities\n", num_cities - missing, num_cities);
  printf("  Total error: %.3f km\n", total_dist);
  printf("  Max error: %.3f km\n", max_dist);
  printf("  Min error: %.3f km\n", min_dist);
  printf("  Average error: %.3f km\n", avg_dist);
  printf("  Median error: %.3f km\n", med_dist);

  delete [] dists;
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
  float* latitudes;
  float* longitudes;

  read_edge(graph_file, num_verts, num_edges, srcs, dsts);
  read_vert_latlong(latlong_graph_file, latitudes, longitudes);
  create_csr(num_verts, num_edges, srcs, dsts, out_array, out_degree_list);
  graph g = {num_verts, num_edges, out_array, out_degree_list, latitudes, longitudes};

  double* CI = centrality_index(&g);
  int* y_hat = match_by_population(&g, CI, populations, num_cities);

  calculate_error(&g, y_hat, latlong, num_cities);

  delete [] srcs;
  delete [] dsts;
  delete [] out_array;
  delete [] out_degree_list;
  delete [] latitudes;
  delete [] longitudes;
  delete [] CI;
  delete [] y_hat;

  delete [] populations;
  for (unsigned int i=0; i<num_cities; ++i) {
    delete [] latlong[i];
  }
  delete [] latlong;
  delete [] cities;

  return EXIT_SUCCESS;
}
