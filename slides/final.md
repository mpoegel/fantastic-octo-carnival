## Final Report

### Transportation Network Analysis
by Matt Poegel


## The Goal

**To use centrality measures such as betweenness and closeness to identify population centers from
road networks and match them with the major cities**


## Data
Roads are undirected edges and road intersections are vertices.

* Minnesota (2,635 nodes): \
  [http://www.cise.ufl.edu/research/sparse/matrices/Gleich/minnesota.html](http://www.cise.ufl.edu/research/sparse/matrices/Gleich/minnesota.html)
* Colorado (435,666 nodes): \
  [http://www.dis.uniroma1.it/challenge9/download.shtml](http://www.dis.uniroma1.it/challenge9/download.shtml)
* Northeast United States (1,524,453 nodes): \
  [http://www.dis.uniroma1.it/challenge9/download.shtml](http://www.dis.uniroma1.it/challenge9/download.shtml)

### Ground Truth Data
* City names plus longitude and latitude: \
  [http://www.realestate3d.com/gps/latlong.htm](http://www.realestate3d.com/gps/latlong.htm)
* City populations for cities in every US state: \
  2015 US Census estimates


## The Big Picture

![Minnesota road network](http://www.cise.ufl.edu/research/sparse/matrices/Gleich/minnesota_gplot_big.png)


## Sample Cities

![Minnesota sample cities](reports/figures/mn_cities.png)


### Ground Truth Data
* City names plus longitude and latitude: \
  [http://www.realestate3d.com/gps/latlong.htm](http://www.realestate3d.com/gps/latlong.htm)
* City populations for cities in every US state: \
  2015 US Census estimates


## The Project

This project is available on my GitHub at [https://github.com/mpoegel/fantastic-octo-carnival](https://github.com/mpoegel/fantastic-octo-carnival).

### Structure
* `./data`
    * `/external`
    * `/processed`
    * `/raw`
* `./setup`- scripts to download datasets and libraries
* `./slides`
* `./src`
    * `graph.h`, `graph.cpp`
    * `utils.h`, `utils.cpp`
    * `main.cpp`
    * `make-datasets.py`
    * `map.R`
* `Makefile`- `make`, `make slides`, `make data`


## The Analysis

### Centrality algorithms on the raw graph 
### Centrality algorithms on the coarsened graph
### Graph kernel methods


## Centrality Algorithms

Ran the analysis with unweighted edges and weighted edges.

* betweenness
    * $C(v_i) = \sum_{j \not= i} \sum_{k \not= i, k > j} \frac{\eta_{jk}(v_i)}{\eta_{jk}}$
    where $\eta_{jk}$ is number of shortest paths between $v_j$ and $v_k$ and \\
          $\eta_{jk}(v_i)$ is the number of such paths that contain $v_i$
* closeness
    * $C(v_i) = \frac{1}{\sum_{j}d(v_i, v_j)}$
* eccentricity
    * $C(v_i) = \frac{1}{\max_j{d(v_i, v_j)}}$
* pagerank


## Results (raw graph)

```
**** analysis on raw graph ****
Betweenness finished after 2.765 seconds
Closeness finished after 0.036 seconds
Pagerank finished after 0.001 seconds
Eccentricity finished after 0.015 seconds
Analysis on 27/28 cities
  Total error: 7678.062 km
  Max error: 568.061 km
  Min error: -0.000 km
  Average error: 274.217 km
  Median error: 295.999 km
```


## Results (raw graph)

![Raw graph analysis: centrality index](reports/figures/mn_raw_centrality_index.png)


## Results (raw graph)

![Raw graph analysis: ground truth vs. estimate](reports/figures/mn_raw_ground_truth_vs_estimates.png)


## Graph Coarsening

Followed the same approach from Homework 3:

* label propagation algorithm
* graph coarsening
* computed average latitudes and longitudes
* used the coarsened weights instead of creating new weights


## Results (coarsened graph)

```
**** analysis on coarsened graph ****
Label Propagation finished after 0.506 seconds
Coarsened graph has 812 edges between 206 communities
  after 0.016 seconds
Betweenness finished after 0.001 seconds
Closeness finished after 0.001 seconds
Pagerank finished after 0.001 seconds
Eccentricity finished after 0.000 seconds
Analysis on 27/28 cities
  Total error: 5266.040 km
  Max error: 401.468 km
  Min error: -0.000 km
  Average error: 188.073 km
  Median error: 210.951 km
```


## Results (raw graph)

![Coarsened graph analysis: centrality index](reports/figures/mn_coarsened_centrality_index.png)


## Results (coarsened graph)

![Coarsened graph analysis: ground truth vs. estimate](reports/figures/mn_coarsened_ground_truth_vs_estimates.png)


## Graph Kernel Methods

### Exponential Diffusion Kernel

Given:

* the adjacency matrix, $\mathbf{A}$
* the diagonal degree matrix, $\Delta$
* the Laplacian matrix, $\mathbf{L} = \mathbf{A} - \Delta$

We choose $\mathbf{S} = \mathbf{L}$ and define the exponential diffusion kernel as follows:

$$\mathbf{K} = \sum^{\infty}_{l=0} \frac{1}{l!}\beta^l\mathbf{S}^l$$


## Graph Kernel Methods

### Kernel K-means
Similar to k-means but uses a kernel

* Create the kernel
    * Expensive operation because this requires the eigendecomposition
    * Used two external libraries: Eigen and Spectra 
* Run Kernel K-means to obtain ~`num_cities` clusters
* Computer the average latitude, longitude, and centrality index for each cluster

## Results (Kernel K-means)

```
**** analysis using kernel-based methods ****
Exponential diffusion kernel computed after 0.385 seconds
K-means (k=28) finished after 3 iterations in 0.035
  seconds
Analysis on 24/28 cities
  Total error: 5042.119 km
  Max error: 449.319 km
  Min error: 0.000 km
  Average error: 180.076 km
  Median error: 183.518 km
```


## Results (Kernel K-means)

![Kernel K-means graph analysis: centrality index](reports/figures/mn_kernel_kmeans_centrality_index.png)


## Results (Kernel K-means)

![Kernel K-means graph analysis: ground truth vs. estimate](reports/figures/mn_kernel_kmeans_ground_truth_vs_estimates.png)


## Conclusion

* Underlying assumpting confirmed
    * Matching graph centrality with population centers is possible
* Matching graph centrality to specific populations centers is very hard
* Coarse graph analysis was better than the raw graph analysis
* Kernel K-means analysis was more conservation in its estimates
* Plenty of room to continue this analysis


## References

* Annual Estimates of the Resident Population: April 1, 2010 to July 1, 2015 Source: U.S. Census Bureau, 
May 2016.
* Guennebaud, G. & Jacob, B. (2010). Eigen (Version 3.3.1).
* Qiu, Y. (2015). Spectra (Version 0.4.0). West Lafayette, IN.
* The University of Florida Sparse Matrix Collection, T. A. Davis and Y. Hu, ACM Transactions on 
Mathematical Software, Vol 38, Issue 1, 2011, pp 1:1 - 1:25. http://www.cise.ufl.edu/research/sparse/matrices
* Zaki, M. & Meira, W. (2014). Data Mining and Analysis: Fundamental Concepts and Algorithms 
(1st ed.). Cambridge University Press.
