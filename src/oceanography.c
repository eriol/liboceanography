#include <math.h>

#include "oceanography.h"

#define A(xt) (-3.107E-3 * (xt) + 0.4215)
#define B(xt) ((4.464e-4 * (xt) + 3.426e-2) * (xt) + 1.0)
#define C(xp) (((3.989e-15 * (xp) - 6.370e-10) * (xp) + 2.070e-5) * (xp))
#define RT35(xt) ((((1.0031e-9 * (xt) - 6.9698e-7) * (xt) + 1.104259e-4) \
                  * (xt) + 2.00564e-2) * (xt) + 0.6766097)

static double _sal(double xr, double xt)
{
    return ((((2.7081 * xr - 7.0261) * xr + 14.0941) * xr + 25.3851)
            * xr - 0.1692) * xr + 0.0080 + (xt / (1.0 + 0.0162 * xt)) *
           (((((-0.0144 * xr + 0.0636) * xr - 0.0375) * xr - 0.0066)
            * xr -0.0056) * xr + 0.0005);
}

static double _dsal(double xr, double xt)
{
    return ((((13.5405 * xr - 28.1044) * xr + 42.2823) * xr + 50.7702)
            * xr - 0.1692) + ( xt / (1.0 + 0.0162 * xt)) *
           ((((-0.0720 * xr + 0.2544) * xr -0.1125) * xr - 0.0132)
            * xr -0.0056);
}

double salinity(double conductivity, double temperature, double pressure)
{
    double corrected_temperature, rt;
    corrected_temperature = temperature - 15.0;

    if (conductivity <= 5e-4)
        return 0.0;

    rt = conductivity / (RT35(temperature) * (1.0 + C(pressure) /
                                              (B(temperature) + A(temperature)
                                              * conductivity)));
    rt = sqrt(fabs(rt));

    return _sal(rt, corrected_temperature);
}

double conductivity(double salinity, double temperature, double pressure)
{
    double corrected_temperature, rt, si, dels, rtt, cp, bt, r;
    int n = 0;
    corrected_temperature = temperature - 15.0;

    if (salinity <= 0.02)
        return 0.0;

    rt = sqrt(salinity / 35.0);
    si = _sal(rt, corrected_temperature);

    do {
        rt = rt + (salinity - si) / _dsal(rt, corrected_temperature);
        si = _sal(rt, corrected_temperature);
        dels = fabs(si - salinity);
        n++;
    } while(n < 10 && dels > 1.0e-4);

    rtt = RT35(temperature) * rt * rt;
    cp = rtt * (C(pressure) + B(temperature));
    bt = B(temperature) - rtt * A(temperature);
    r = sqrt(fabs(bt * bt + 4.0 * A(temperature) * cp)) - bt;

    return 0.5 * r / A(temperature);
}
