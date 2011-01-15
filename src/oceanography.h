/*
 * Copyright 2011 Daniele Tricoli <eriol@mornie.org>
 *
 * This file is part of liboceanography.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 3 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library. If not, see http://www.gnu.org/licenses/.
 *
 *
 * oceanography.h -- liboceanography header.
 */

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

#endif /* OCEANOGRAPHY_H */
