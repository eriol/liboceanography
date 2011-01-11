#include <math.h>
#include <stdlib.h>

#include <check.h>

#include "../src/oceanography.h"

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

Suite *oceanography_suite(void)
{
    Suite *s = suite_create("Oceanography");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_salinity);
    tcase_add_test(tc_core, test_conductivity);
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
