#ifndef OCEANOGRAPHY_H
#define OCEANOGRAPHY_H

double salinity(double conductivity, double temperature, double pressure);
double conductivity(double salinity, double temperature, double pressure);
double specific_volume_anomaly(double salinity, double temperature,
                               double pressure, double *sigma);
double depth(double pressure, double latitude);
double freezing_point(double salinity, double pressure);
double specific_heat(double salinity, double temperature, double pressure);
double adiabatic_temperature_gradient(double salinity, double temperature,
                                      double pressure);
double potential_temperature(double salinity, double temperature,
                             double pressure, double reference_pressure);
double sound_speed(double salinity, double temperature, double pressure);

#endif
