#! /bin/bash

mkdir -p lib/Eigen
wget http://bitbucket.org/eigen/eigen/get/3.3.1.tar.bz2 -P lib
tar -vxjf lib\3.3.1.tar.bz2 -C lib\Eigen --strip-components=1
rm -f lib\3.3.1.tar.bz2
