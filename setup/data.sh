#! /bin/bash

wget http://www.cise.ufl.edu/research/sparse/MM/Gleich/minnesota.tar.gz -P data/raw
tar -zxvf data/raw/minnesota.tar.gz -C data/raw
rm -rf data/raw/minnesota.tar.gz

wget http://www.realestate3d.com/gps/latlong.htm -O data/raw/latlong-cities.htm
