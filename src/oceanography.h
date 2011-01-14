#ifndef OCEANOGRAPHY_H
#define OCEANOGRAPHY_H

double salinity(double conductivity, double temperature, double pressure);
double conductivity(double salinity, double temperature, double pressure);
double specific_volume_anomaly(double salinity, double temperature,
                               double pressure, double *sigma);

#endif
