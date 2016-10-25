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
