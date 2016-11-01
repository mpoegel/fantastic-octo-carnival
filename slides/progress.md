## Progress Report

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


## Progress
* Downloaded all of the data for the Minnesota network
* Parse all five input data files: two graph files and three ground truth files
* Calculate the betweenness centrality measure
* Normalize the centrality index and populations to be between 0 and 1
* For each city, find the vertex with the closest centrality index: $\mathbf{\hat{y}}$
* Calculate the distance between $\mathbf{\hat{y}}$ and $\mathbf{y}$ in kilometers
    * report the total, max, min, average, and median error


## Intermediate Results
* High error so far with Minnesota network, not sure why
```
Betweenness finished after 70.355 seconds
Analysis on 27/28 cities
  Total error: 5924.000 km
  Max error: 482.000 km
  Min error: -0.000 km
  Average error: 211.571 km
  Median error: 195.500 km
```


## Next Steps
* Determine why Minnesota's results are not very good
    * Consider plotting $\mathbf{\hat{y}}$ against $\mathbf{y}$ and distribution of errors
* Parallelize everything that can be
* Add more centrality measures and combine them in some way to form a centrality index
* Create the ground truth data sets for Colorado and the Northeast and run algorithm
