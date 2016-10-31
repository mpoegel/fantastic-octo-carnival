#include "graph.h"


void read_edge(char* filename, int &num_verts, int &num_edges, int* &srcs, int* &dsts)
{
  ifstream infile;
  string line;
  infile.open(filename);

  getline(infile, line, '\n');
  while (line[0] == '%') {
    getline(infile, line, '\n');
  }
  istringstream iss(line);
  string tmp;
  iss >> tmp;
  num_verts = atoi(tmp.c_str());
  iss >> tmp;
  num_edges = atoi(tmp.c_str());
  iss >> tmp;

  int src, dst;
  int counter = 0;

  num_edges *= 2;
  srcs = new int[num_edges];
  dsts = new int[num_edges];
  for (unsigned int i=0; i<num_edges / 2; i++) {
    getline(infile, line, ' ');
    src = atoi(line.c_str());
    getline(infile, line);
    dst = atoi(line.c_str());

    srcs[counter] = src;
    dsts[counter] = dst;
    counter++;
    srcs[counter] = dst;
    dsts[counter] = src;
    counter++;
  }

  infile.close();
}


void read_vert_latlong(char* filename, float* &latitudes, float* &longitudes)
{
  ifstream infile;
  string line;
  infile.open(filename);

  getline(infile, line, '\n');
  while (line[0] == '%') {
    getline(infile, line, '\n');
  }

  istringstream iss(line);
  string tmp;
  iss >> tmp;
  int num_verts = atoi(tmp.c_str());
  iss >> tmp;

  latitudes = new float[num_verts];
  longitudes = new float[num_verts];

  for (int i=0; i<num_verts; ++i) {
    getline(infile, line);
    longitudes[i] = atof(line.c_str());
  }
  for (int i=0; i<num_verts; ++i) {
    getline(infile, line);
    latitudes[i] = atof(line.c_str());
  }

  infile.close();
}


void create_csr(int num_verts, int num_edges, int* srcs, int* dsts, int* &out_array,
                int* &out_degree_list)
{
  out_array = new int[num_edges];
  out_degree_list = new int[num_verts+1];

  for (unsigned int i=0; i<num_edges; i++) {
    out_array[i] = 0;
  }
  for (unsigned int i=0; i<num_verts + 1; i++) {
    out_degree_list[i] = 0;
  }

  int* temp_counts = new int[num_verts];
  for (unsigned int i=0; i<num_verts; i++) {
    temp_counts[i] = 0;
  }
  for (unsigned int i=0; i<num_edges; i++) {
    temp_counts[srcs[i]]++;
  }
  for (unsigned int i=0; i<num_verts; i++) {
    out_degree_list[i + 1] = out_degree_list[i] + temp_counts[i];
  }
  copy(out_degree_list, out_degree_list + num_verts, temp_counts);
  for (unsigned int i=0; i<num_edges; i++) {
    out_array[temp_counts[srcs[i]]++] = dsts[i];
  }
  
  delete [] temp_counts;
}


double* centrality_index(graph* g)
{
  // calculate individual measures
  int* betweenness = calc_betweenness(g);

  // aggregate measures together to form the final index and normalize
  double* CI = new double[g->num_verts];
  double max = 0.0;
  for (unsigned int i=0; i<g->num_verts; ++i) {
    double val = (double)betweenness[i];
    CI[i] = val;
    if (val > max) {
      max = val;
    }
  }
  for (unsigned int i=0; i<g->num_verts; ++i) {
    CI[i] /= max;
  }

  delete [] betweenness;

  return CI;
}


int* calc_betweenness(graph* g)
{
  double timer = omp_get_wtime();

  // get the shortest paths using Floyd-Warshall
  int** dist = new int*[g->num_verts];
  int** next = new int*[g->num_verts];
  for (int i=0; i<g->num_verts; ++i) {
    dist[i] = new int[g->num_verts];
    next[i] = new int[g->num_verts];
    for (int k=0; k<g->num_verts; ++k) {
      dist[i][k] = -1;
      next[i][k] = -1;
    }
    dist[i][i] = 0;
    next[i][i] = i;
  }

  for (int v=0; v<g->num_verts; ++v) {
    int out_degree = out_degree(g, v);
    int* out_vertices = out_vertices(g, v);
    for (int k=0; k<out_degree; ++k) {
      int u = out_vertices[k];
      dist[v][u] = 1;
      next[v][u] = u;
    }
  }

  for (int k=0; k<g->num_verts; ++k) {
    for (int i=0; i<g->num_verts; ++i) {
      for (int j=0; j<g->num_verts; ++j) {
        if (dist[i][k] == -1 || dist[k][j] == -1) {
          continue;
        }
        if (dist[i][j] == -1 || dist[i][k] + dist[k][j] < dist[i][j]) {
          dist[i][j] = dist[i][k] + dist[k][j];
          next[i][j] = next[i][k];
        }
      }
    }
  }
  
  int* betweenness = new int[g->num_verts];
  // count the shortest paths that go through each vertex
  for (int i=0; i<g->num_verts; ++i) {
    betweenness[i] = 0;
  }
  for (int u=0; u<g->num_verts; ++u) {
    for (int v=u+1; v<g->num_verts; ++v) {
      int s = u;
      if (dist[s][v] == -1) {
        continue;
      }
      while (s != v) {
        s = next[s][v];
        betweenness[s]++;
      }
    }
  }

  timer = omp_get_wtime() - timer;
  printf("Betweenness finished after %.3f seconds\n", timer);

  for (int i=0; i<g->num_verts; ++i) {
    delete [] dist[i];
    delete [] next[i];
  }
  delete [] dist;
  delete [] next;

  return betweenness;
}


bool betweenness_comp(int* a, int* b)
{
  return a[1] > b[1];
}


void output_info(graph* g, int* betweenness)
{
  int** centrality_map = new int*[g->num_verts];
  for (int i=0; i<g->num_verts; ++i) {
    centrality_map[i] = new int[2];
    centrality_map[i][0] = i;
    centrality_map[i][1] = betweenness[i];
  }
  printf("Highest betweenness:\n");
  sort(centrality_map, centrality_map + g->num_verts, betweenness_comp);
  for (int i=0; i<5; ++i) {
    printf("  %d - %d\n", centrality_map[i][0], centrality_map[i][1]);
  }
}


int* match_by_population(graph* g, double* CI, int* populations, int num_cities)
{
  // normalize the populations
  double* norm_pops = new double[num_cities];
  unsigned int max = 0;
  for (unsigned int i=0; i<num_cities; ++i) {
    if (populations[i] > 0 && populations[i] > max) {
      max = populations[i];
    }
  }

  int* y_hat = new int[num_cities];
  for (unsigned int i=0; i<num_cities; ++i) {
    norm_pops[i] = (double)populations[i] / (double)max;
    double min_dist = INT_MAX;
    double d;
    unsigned int min_v = -1;
    if (norm_pops[i] > 0) {
      for (unsigned int v=0; v<g->num_verts; ++v) {
        d = fabs(norm_pops[i] - CI[v]);
        if (d < min_dist) {
          min_dist = d;
          min_v = v;
        }
      }
    }
    y_hat[i] = min_v;
  }

  delete [] norm_pops;

  return y_hat;
}