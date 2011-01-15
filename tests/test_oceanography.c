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
 * test_oceanography.c -- Unit tests for liboceanography.
 */

#include <math.h>
#include <stdlib.h>

#include <check.h>

#include "oceanography.h"

#define EPSILON 0.00001


static int cmp_double(double x, double y)
{
    return fabs(x - y) < EPSILON;
}


START_TEST(test_salinity)
{
    ck_assert(cmp_double(salinity(1, 15, 0), 35.0));
    ck_assert(cmp_double(salinity(1.2, 20, 2000), 37.245628));
    ck_assert(cmp_double(salinity(0.65, 5, 1500), 27.995347));
    ck_assert(cmp_double(salinity(1.888091, 40, 10000), 40.0));

    ck_assert(cmp_double(salinity(5e-5, 15, 0), 0.0));
}
END_TEST

START_TEST(test_conductivity)
{
    ck_assert(cmp_double(conductivity(35.0, 15, 0), 1.0));
    ck_assert(cmp_double(conductivity(37.245628, 20, 2000), 1.2));
    ck_assert(cmp_double(conductivity(27.995347, 5, 1500), 0.65));
    ck_assert(cmp_double(conductivity(40.0, 40, 10000), 1.888091));

    ck_assert(cmp_double(conductivity(0.02, 15, 0), 0.0));
}
END_TEST

START_TEST(test_specific_volume_anomaly)
{
    double sigma;

    ck_assert(cmp_double(specific_volume_anomaly(0, 0, 0, &sigma),
                         2749.539368));
    ck_assert(cmp_double(sigma, -0.1574));

    ck_assert(cmp_double(specific_volume_anomaly(0, 0, 1000, &sigma),
                             2692.644915));
    ck_assert(cmp_double(sigma, 4.872729));

    ck_assert(cmp_double(specific_volume_anomaly(40, 0, 0, &sigma),
                         -380.789102));
    ck_assert(cmp_double(sigma, 32.147101));

    ck_assert(cmp_double(specific_volume_anomaly(40, 40, 10000, &sigma),
                         981.301907));
    ck_assert(cmp_double(sigma, 59.820375));
}
END_TEST

START_TEST(test_depth)
{
    ck_assert(cmp_double(depth(500, 0), 496.652992));
    ck_assert(cmp_double(depth(10000, 30), 9712.653072));
    ck_assert(cmp_double(depth(10000, 90), 9674.231441));
}
END_TEST

START_TEST(test_freezing_point)
{
    ck_assert(cmp_double(freezing_point(5, 0), -0.273763));
    ck_assert(cmp_double(freezing_point(20, 300), -1.309106));
    ck_assert(cmp_double(freezing_point(40, 500), -2.588567));
}
END_TEST

START_TEST(test_specific_heat)
{
    ck_assert(cmp_double(specific_heat(25, 0, 0), 4048.440412));
    ck_assert(cmp_double(specific_heat(35, 20, 5000), 3894.992770));
    ck_assert(cmp_double(specific_heat(40, 40, 10000), 3849.499481));
}
END_TEST

START_TEST(test_adiabatic_temperature_gradient)
{
    ck_assert(cmp_double(adiabatic_temperature_gradient(25, 0, 0),
                         1.687100e-05));
    ck_assert(cmp_double(adiabatic_temperature_gradient(30, 20, 9000),
                         2.416120e-04));
    ck_assert(cmp_double(adiabatic_temperature_gradient(40, 40, 10000),
                         3.255976e-04));
}
END_TEST

START_TEST(test_potential_temperature)
{
    ck_assert(cmp_double(potential_temperature(25, 0, 0, 0), 0));
    ck_assert(cmp_double(potential_temperature(25, 40, 0, 0), 40.000000));
    ck_assert(cmp_double(potential_temperature(30, 20, 9000, 0),
                         18.296537));
    ck_assert(cmp_double(potential_temperature(40, 40, 10000, 0),
                         36.998992));
}
END_TEST

START_TEST(test_sound_speed)
{
    ck_assert(cmp_double(sound_speed(25, 0, 0), 1435.789875));
    ck_assert(cmp_double(sound_speed(35, 20, 5000), 1604.476282));
    ck_assert(cmp_double(sound_speed(40, 40, 10000), 1731.995394));
}
END_TEST


Suite *oceanography_suite(void)
{
    Suite *s = suite_create("Oceanography");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_salinity);
    tcase_add_test(tc_core, test_conductivity);
    tcase_add_test(tc_core, test_specific_volume_anomaly);
    tcase_add_test(tc_core, test_depth);
    tcase_add_test(tc_core, test_freezing_point);
    tcase_add_test(tc_core, test_specific_heat);
    tcase_add_test(tc_core, test_adiabatic_temperature_gradient);
    tcase_add_test(tc_core, test_potential_temperature);
    tcase_add_test(tc_core, test_sound_speed);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = oceanography_suite();
  SRunner *sr = srunner_create(s);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
