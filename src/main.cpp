#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

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


int* random_assignment(graph* g, int num_cities)
{
  int* y_hat = new int[num_cities];
  for (unsigned int i=0; i<num_cities; ++i) {
    y_hat[i] = rand() % g->num_verts;
  }
  return y_hat;
}


void calculate_error(float* g_lat, float* g_lon, int* y_hat, float** latlong, int num_cities)
{
  double* dists = new double[num_cities];
  double total_dist = 0.0;
  unsigned int missing = 0;
  for (unsigned int i=0; i<num_cities; ++i) {
    int v = y_hat[i];
    if (v > 0) {
      double d = measureLatLongDist(g_lat[v], g_lon[v], latlong[i][0],
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


void save_yhat(const char* filename, float* g_lats, float* g_lons, int* y_hat, int num_cities)
{
  FILE* fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: could not open output file %s\n", filename);
    return;
  }
  fprintf(fp, "id,lat,lon\n");
  float epsilon = 0.001;
  for (unsigned int i=0; i<num_cities; ++i) {
    int v = y_hat[i];
    if (abs(g_lats[v]) < epsilon || abs(g_lons[v]) < epsilon) {
      continue;
    }
    fprintf(fp, "%d,%f,%f\n", i, g_lats[v], g_lons[v]);
  }
  fclose(fp);
}


void save_CI(const char* filename, float* g_lats, float* g_lons, unsigned int N, double* CI)
{
  FILE* fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: could not open output file %s\n", filename);
    return;
  }
  fprintf(fp, "id,lat,lon,CI\n");
  for (unsigned int v=0; v<N; ++v) {
    if (CI[v] > 0) {
      fprintf(fp, "%d,%f,%f,%f\n", v, g_lats[v], g_lons[v], CI[v]);
    } 
  }
  fclose(fp);
}

void save_clusters(const char* filename, graph* g, int* clusters)
{
  FILE* fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: could not open output file %s\n", filename);
    return;
  }
  fprintf(fp, "id,lat,lon,cluster\n");
  for (unsigned int v=0; v<g->num_verts; ++v) {
    fprintf(fp, "%d,%f,%f,%d\n", v, g->latitudes[v], g->longitudes[v], clusters[v]);
  }
  fclose(fp);
}


int main(int argc, char* argv[])
{

  if (argc < 6) {
    fprintf(stderr, "Usage: %s <input_graph> <latlong_graph> <pop_file> <latlong_file>"
            "<cities_file> [output_file]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  srand(time(0));

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
  double* out_weights;
  int* out_degree_list;
  float* latitudes;
  float* longitudes;
  int max_degree;

  read_edge(graph_file, num_verts, num_edges, srcs, dsts);
  read_vert_latlong(latlong_graph_file, latitudes, longitudes);
  create_csr(num_verts, num_edges, srcs, dsts, latitudes, longitudes, out_array, out_weights,
             out_degree_list, max_degree);
  graph g = {num_verts, num_edges, out_array, out_degree_list, latitudes, longitudes, out_weights,
             max_degree};
  double* CI;
  int* y_hat;

  /*
   *
   */
  printf("**** random assignment baseline ****\n");

  y_hat = random_assignment(&g, num_cities);

  calculate_error(g.latitudes, g.longitudes, y_hat, latlong, num_cities);

  delete[] y_hat;


  /*
   * Run analysis on the graph as it is
   */
  printf("**** analysis on raw graph ****\n");

  CI = centrality_index(&g);
  y_hat = match_by_population(CI, g.num_verts, populations, num_cities);

  calculate_error(g.latitudes, g.longitudes, y_hat, latlong, num_cities);

  if (argc >= 7) {
    string base_filename = string(argv[6]);
    const char* yhat_filename = string(base_filename + ".yhat").c_str();
    save_yhat(yhat_filename, g.latitudes, g.longitudes, y_hat, num_cities);
    const char* ci_filename = string(base_filename + ".ci").c_str();    
    save_CI(ci_filename, g.latitudes, g.longitudes, g.num_verts, CI);
  }

  delete[] CI;
  delete[] y_hat;

  /*
   * Run analysis on coarsened graph
   */
  printf("\n**** analysis on coarsened graph ****\n");
  
  int* labels;
  graph coarse_g;
  coarsen(&g, &coarse_g, labels);
  
  CI = centrality_index(&coarse_g);

  y_hat = match_by_population(CI, g.num_verts, populations, num_cities);

  calculate_error(coarse_g.latitudes, coarse_g.longitudes, y_hat, latlong, num_cities);

  if (argc >= 7) {
    string base_filename = string(argv[6]);
    const char* yhat_filename = string(base_filename + ".coarse.yhat").c_str();
    save_yhat(yhat_filename, coarse_g.latitudes, coarse_g.longitudes, y_hat, num_cities);
    const char* ci_filename = string(base_filename + ".coarse.ci").c_str();    
    save_CI(ci_filename, coarse_g.latitudes, coarse_g.longitudes, coarse_g.num_verts, CI);
  }

  delete [] labels;
  delete [] CI;
  delete [] y_hat;


  /*
   * Graph kernel methods
   */
  printf("\n**** analysis using kernel-based methods ****\n");   

  MatrixXd K(g.num_verts, g.num_verts);
  double beta = 1.0; 
  exponential_diffusion_kernel(&coarse_g, K, beta);
  unsigned int num_clusters = num_cities;
  int* clusters = kernel_kmeans(K, num_clusters, 0.1);

  double* cluster_CI = average_by_cluster(clusters, num_clusters, CI, coarse_g.num_verts);
  double max_ci = 0;
  for (unsigned int i=0; i<num_clusters; ++i) if (cluster_CI[i] > max_ci) max_ci = cluster_CI[i];
  for (unsigned int i=0; i<num_clusters; ++i) cluster_CI[i] /= max_ci;
  y_hat = match_by_population(cluster_CI, num_clusters, populations, num_cities);
  float* cluster_lat = average_by_cluster(clusters, num_clusters, coarse_g.latitudes,
                                          coarse_g.num_verts);
  float* cluster_lon = average_by_cluster(clusters, num_clusters, coarse_g.longitudes,
                                          coarse_g.num_verts);
  calculate_error(cluster_lat, cluster_lon, y_hat, latlong, num_cities);

  if (argc >= 7) {
    string base_filename = string(argv[6]);
    const char* yhat_filename = string(base_filename + ".coarse.clusters.yhat").c_str();
    save_yhat(yhat_filename, cluster_lat, cluster_lon, y_hat, num_cities);
    const char* cluster_filename = string(base_filename + ".coarse.clusters").c_str();
    save_clusters(cluster_filename, &coarse_g, clusters);
    const char* ci_filename = string(base_filename + ".coarse.clusters.ci").c_str();
    save_CI(ci_filename, cluster_lat, cluster_lon, num_clusters, cluster_CI);
  }

  delete[] clusters;

  /*
   * Cleanup
   */

  delete [] srcs;
  delete [] dsts;
  delete [] out_array;
  delete [] out_degree_list;
  delete [] latitudes;
  delete [] longitudes;

  delete [] populations;
  for (unsigned int i=0; i<num_cities; ++i) {
    delete [] latlong[i];
  }
  delete [] latlong;
  delete [] cities;

  return EXIT_SUCCESS;
}
