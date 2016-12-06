#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

#define RADIUS_EARTH 6378.137


double measureLatLongDist(double lat1, double lon1, double lat2, double lon2);


template<typename T>
void displayn(T* arr, unsigned int n)
{
  for (unsigned int i=0; i<n; ++i) {
    cout << i << ": " << arr[i] << endl; 
  }
}


template<typename T>
void fdisplayn(T* arr, unsigned int n, unsigned int sig)
{
  cout << setprecision(sig);
  for (unsigned int i=0; i<n; ++i) {
    cout << i << ": " << arr[i] << endl; 
  }
}


template<typename T>
void displaynd(T* arr, unsigned int n, unsigned int d)
{
  for (unsigned int i=0; i<n; ++i) {
    cout << i << ": "; 
    for (unsigned int k=0; k<d; ++k) {
      cout << arr[i][k] << " ";
    }
    cout << endl;
  }
}


template<typename T>
void fdisplaynd(T arr, unsigned int n, unsigned int d, unsigned int sig)
{
  cout << setprecision(sig);
  for (unsigned int i=0; i<n; ++i) {
    cout << i << ": "; 
    for (unsigned int k=0; k<d; ++k) {
      cout << arr[i][k] << " ";
    }
    cout << endl;
  }
}

#endif
