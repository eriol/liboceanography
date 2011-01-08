#ifndef OCEANOGRAPHY_H
#define OCEANOGRAPHY_H

double salinity(double conductivity, double temperature, double pressure);
double conductivity(double salinity, double temperature, double pressure);

#endif