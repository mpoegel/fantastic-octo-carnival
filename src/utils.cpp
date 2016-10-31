#include "utils.h"


int measureLatLongDist(double lat1, double lon1, double lat2, double lon2)
{
  double dLat = lat2 * M_PI / 180 - lat1 * M_PI / 180;
  double dLon = lon2 * M_PI / 180 - lon1 * M_PI / 180;

  double a = sin(dLat/2) * sin(dLat/2) + cos(lat1 * M_PI / 180) * cos(lat2 * M_PI / 180) *
             sin(dLon/2) * sin(dLon/2);
  double c = 2 * atan2(sqrt(a), sqrt(1-a));
  double d = RADIUS_EARTH * c;
  return d;
}
