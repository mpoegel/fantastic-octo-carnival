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


void create_csr(int num_verts, int num_edges, int* srcs, int* dsts, float* latitudes, 
                float* longitudes, int* &out_array, double* &out_weights, int* &out_degree_list)
{
  out_array = new int[num_edges];
  out_weights = new double[num_edges];
  out_degree_list = new int[num_verts+1];

  for (unsigned int i=0; i<num_edges; i++) {
    out_array[i] = 0;
    out_weights[i] = 0;
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
    out_array[temp_counts[srcs[i]]] = dsts[i];
    double w = measureLatLongDist(latitudes[srcs[i]], longitudes[srcs[i]], latitudes[dsts[i]],
                                  longitudes[dsts[i]]);
    out_weights[temp_counts[srcs[i]]++] = w;
  }
  
  delete [] temp_counts;
}


void shortest_paths(graph* g, double** dist, int** next)
{
  // get the shortest paths using Floyd-Warshall
  for (int i=0; i<g->num_verts; ++i) {
    dist[i] = new double[g->num_verts];
    if (next) next[i] = new int[g->num_verts];
    for (int k=0; k<g->num_verts; ++k) {
      dist[i][k] = -1;
      if (next) next[i][k] = -1;
    }
    dist[i][i] = 0;
    if (next) next[i][i] = i;
  }

  for (int v=0; v<g->num_verts; ++v) {
    int out_degree = out_degree(g, v);
    int* out_vertices = out_vertices(g, v);
    double* out_weights = out_weights(g, v);
    for (int k=0; k<out_degree; ++k) {
      int u = out_vertices[k];
      dist[v][u] = out_weights[k];
      if (next) next[v][u] = u;
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
          if (next) next[i][j] = next[i][k];
        }
      }
    }
  }
}


double* centrality_index(graph* g)
{
  double** dist = new double*[g->num_verts];
  int** next = new int*[g->num_verts];
  shortest_paths(g, dist, next);
  // calculate individual measures
  double* betweenness = calc_betweenness(g, dist, next);
  double* closeness = calc_closeness(g, dist);
  double* pagerank = calc_pagerank(g, 20, 0.85);
  double* eccentricity = calc_eccentricity(g, dist);

  // aggregate measures together to form the final index and normalize
  double* CI = new double[g->num_verts];
  double max = 0.0;
  for (unsigned int i=0; i<g->num_verts; ++i) {
    double b = betweenness[i];
    double c = closeness[i];
    double p = pagerank[i];
    double e = eccentricity[i];
    double val = b; // (0.5 * b) + (0.5 * e);
    CI[i] = val;
    if (isnan(val)) {
      printf("nope nope nope!\n");
    }
    if (val > max) {
      max = val;
    }
  }
  for (unsigned int i=0; i<g->num_verts; ++i) {
    CI[i] /= max;
  }

  delete [] betweenness;
  delete [] closeness;
  delete [] pagerank;
  delete [] eccentricity;
  
  for (int i=0; i<g->num_verts; ++i) {
    delete [] dist[i];
    delete [] next[i];
  }
  delete [] dist;
  delete [] next;

  return CI;
}


double* calc_betweenness(graph* g, double** dist, int** next)
{
  double timer = omp_get_wtime();
  
  double* betweenness = new double[g->num_verts];
  // count the shortest paths that go through each vertex
  for (int i=0; i<g->num_verts; ++i) {
    betweenness[i] = 0.0;
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

  double max_b = 0;
  for (int i=0; i<g->num_verts; ++i) {
    if (betweenness[i] > max_b) {
      max_b = betweenness[i];
    }
  }
  for (int i=0; i<g->num_verts; ++i) {
    betweenness[i] /= max_b;
  }

  timer = omp_get_wtime() - timer;
  printf("Betweenness finished after %.3f seconds\n", timer);

  return betweenness;
}


double* calc_closeness(graph* g, double** dist)
{
  double timer = omp_get_wtime();

  double* closeness = new double[g->num_verts];

  for (int i=0; i<g->num_verts; ++i) {
    closeness[i] = 0.0;
    for (int k=0; k<g->num_verts; ++k) {
      if (dist[i][k] > 0) {
        closeness[i] += dist[i][k];
      }
    }
    closeness[i] = 1.0 / closeness[i];
  }

  timer = omp_get_wtime() - timer;
  printf("Closeness finished after %.3f seconds\n", timer);

  return closeness;
}


#define DELTA 1E-7

double* calc_pagerank(graph* g, unsigned int num_iter, double damping_factor)
{
  double* pageranks = new double[g->num_verts];
  double* pagerank_next = new double[g->num_verts];
  double sum_sinks = 0.0;
  double sum_sinks_next = 0.0;
  
  double timer = omp_get_wtime();
  for (int vert = 0; vert < g->num_verts; ++vert) {
    pageranks[vert] = 1 / (double)g->num_verts;
    int out_degree = out_degree(g, vert);
    if (out_degree > 0)
      pageranks[vert] /= (double)out_degree;
    else
    {
      pageranks[vert] /= (double)g->num_verts;
      sum_sinks += pageranks[vert];
    }
  }

  int* work_queue = new int[g->num_verts];
  int* queue_next = new int[g->num_verts];
  int queue_size = g->num_verts;
  int next_size = 0;

  int num_updates = 0;
  #pragma omp parallel
  {
  #pragma omp for
  for (int vert = 0; vert < g->num_verts; ++vert)
    work_queue[vert] = vert;

  for (int iter = 0; iter < num_iter; ++iter) {
    #pragma omp for reduction(+:num_updates)
    for (int i = 0; i < queue_size; ++i) {
      int vert = work_queue[i];
      double new_pagerank = sum_sinks / (double)g->num_verts;

      int out_degree = out_degree(g, vert);
      int* out_vertices = out_vertices(g, vert);
      for (int j = 0; j < out_degree; ++j) 
        new_pagerank += pageranks[out_vertices[j]];

      new_pagerank *= damping_factor;
      new_pagerank += ((1.0 - damping_factor) / (double)g->num_verts);
      out_degree = out_degree(g, vert);
      if (out_degree > 0)
        new_pagerank /= (double) out_degree;
      else {
        new_pagerank /= g->num_verts;
        sum_sinks_next += pageranks[vert];
      }

      if (fabs(new_pagerank - pageranks[vert]) < DELTA ) {
        ++num_updates;
        int index = 0;
        #pragma omp atomic capture
        index = next_size++;
        queue_next[index] = vert;
      }
      pagerank_next[vert] = new_pagerank;
    }

    #pragma omp single
    {
    double* temp = pageranks;
    pageranks = pagerank_next;
    pagerank_next = temp;
    sum_sinks = sum_sinks_next;
    sum_sinks_next = 0.0;

    int* temp1 = work_queue;
    work_queue = queue_next;
    queue_next = temp1;
    queue_size = next_size;
    next_size = 0;
    num_updates = 0;
    }
  } // end iter for
  } // parallel

  for (int vert = 0; vert < g->num_verts; ++vert)
  {
    if (out_degree(g, vert) > 0)
      pageranks[vert] *= (double)out_degree(g, vert);
    else
      pageranks[vert] *= (double)g->num_verts;
  }
  timer = omp_get_wtime() - timer;
  printf("Pagerank finished after %.3f seconds\n", timer);

  delete [] pagerank_next;

  return pageranks;
}


double* calc_eccentricity(graph* g, double** dist)
{
  double timer = omp_get_wtime();

  double* eccentricity = new double[g->num_verts];

  for (int i=0; i<g->num_verts; ++i) {
    double max_dist = 0.0;
    for (int k=0; k<g->num_verts; ++k) {
      if (dist[i][k] > max_dist) {
        max_dist = dist[i][k];
      }
    }
    if (max_dist > 0) {
      eccentricity[i] = 1.0 / max_dist;
    } else {
      eccentricity[i] = 0.0;
    }
  }

  double max_e;
  for (int i=0; i<g->num_verts; ++i) {
    if (eccentricity[i] > max_e) {
      max_e = eccentricity[i];
    }
  }
  for (int i=0; i<g->num_verts; ++i) {
    eccentricity[i] /= max_e;
  }

  timer = omp_get_wtime() - timer;
  printf("Eccentricity finished after %.3f seconds\n", timer);

  return eccentricity;
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
