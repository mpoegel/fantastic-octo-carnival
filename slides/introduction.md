## Introduction

### Transportation Network Analysis
by Matt Poegel


## The Goal

**use centrality measures such as betweenness and closeness to identify population centers from
road networks and match them with the major cities**


### Data
Roads are undirected edges and road intersections are vertices.

* Minnesota (2,635 nodes): \
  [http://www.cise.ufl.edu/research/sparse/matrices/Gleich/minnesota.html](http://www.cise.ufl.edu/research/sparse/matrices/Gleich/minnesota.html)
* Colorado (435,666 nodes): \
  [http://www.dis.uniroma1.it/challenge9/download.shtml](http://www.dis.uniroma1.it/challenge9/download.shtml)
* Northeast United States (1,524,453 nodes): \
  [http://www.dis.uniroma1.it/challenge9/download.shtml](http://www.dis.uniroma1.it/challenge9/download.shtml)

## Ground Truth Data
* City names plus longitude and latitude: \
  [http://www.realestate3d.com/gps/latlong.htm](http://www.realestate3d.com/gps/latlong.htm)
* City populations for cities in every US state: \
  2015 US Census estimates

## The Big Picture

![Minnesota road network](http://www.cise.ufl.edu/research/sparse/matrices/Gleich/minnesota_gplot_big.png)


## Past Work

* Leskovec et al. used a similar data set of California along with many others in a large study of
community structure in large networks. \
[http://arxiv.org/abs/0810.1355](http://arxiv.org/abs/0810.1355)
* Michael T. Chang looked at urban street patterns in a study of the road networks of San
Francisco, CA and Oldenburg, Germany. \
[http://humnet.scripts.mit.edu/wordpress2/wp-content/uploads/2011/09/Chang_1204_report1.pdf](http://humnet.scripts.mit.edu/wordpress2/wp-content/uploads/2011/09/Chang_1204_report1.pdf)


## Expected Outcome
* Create a method for taking an unlabeled network and attaching semantic data to the nodes and
edges based on the graph’s structure
* Results will be evaluated using the longitude and latitude data that accompanies each of the road
networks
    * Average absolute distance error metric
    * Logarithmic loss


## Progress

* Set up a project repository
* Located sources for population and location data 
* Wrote scripts to download data sets and build the ground truth file for Minnesota
