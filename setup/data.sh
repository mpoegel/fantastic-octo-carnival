#! /bin/bash

wget http://www.cise.ufl.edu/research/sparse/MM/Gleich/minnesota.tar.gz -P data/raw
tar -zxvf data/raw/minnesota.tar.gz -C data/raw
rm -rf data/raw/minnesota.tar.gz

wget http://www.realestate3d.com/gps/latlong.htm -O data/raw/latlong-cities.htm
wget http://www.togetherweteach.com/TWTIC/uscityinfo/23mn/23mn.htm -O data/raw/MN-cities-pop.htm
