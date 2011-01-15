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

double specific_volume_anomaly(double salinity, double temperature,
                               double pressure, double *sigma)
{
    double sig, sr, r1, r2, r3, r4, a, b, c, d, e, a1, b1, aw, bw;
    double ko, kw, k35, v350p, sva, gam, pk, dr35p, dk, dvan;

    double r3500 = 1028.1063;
    double dr350 = 28.106331;
    r4 = 4.8314e-4;

    pressure = pressure / 10.;
    sr = sqrt(fabs(salinity));

    r1 = ((((6.536332e-9 * temperature - 1.120083e-6) * temperature +
          1.001685e-4) * temperature - 9.095290e-3) * temperature +
          6.793952e-2) * temperature - 28.263737;
    r2 = (((5.3875e-9 * temperature - 8.2467e-7) * temperature + 7.6438e-5) *
          temperature - 4.0899e-3) * temperature + 8.24493e-1;
    r3 = (-1.6546e-6 * temperature + 1.0227e-4) * temperature - 5.72466e-3;
    sig = (r4 * salinity + r3 * sr + r2) * salinity + r1;
    v350p = 1.0 / r3500;
    sva = -sig * v350p / (r3500 + sig);
    *sigma = sig + dr350;

    if (pressure == 0.0)
        return sva * 1.0e+8;

    e = (9.1697e-10 * temperature + 2.0816e-8) * temperature - 9.9348e-7;
    bw = (5.2787e-8 * temperature - 6.12293e-6) * temperature + 3.47718e-5;
    b = bw + e * salinity;

    d = 1.91075e-4;
    c = (-1.6078e-6 * temperature - 1.0981e-5) * temperature + 2.2838e-3;
    aw = ((-5.77905e-7 * temperature + 1.16092e-4) * temperature +
          1.43713e-3) * temperature - 0.1194975;
    a = (d * sr + c) * salinity + aw;

    b1 = (-5.3009e-4 * temperature + 1.6483e-2) * temperature + 7.944e-2;
    a1 = ((-6.1670e-5 * temperature + 1.09987e-2) * temperature - 0.603459) *
         temperature + 54.6746;
    kw = (((-5.155288e-5 * temperature + 1.360477e-2) * temperature -
          2.327105) * temperature +148.4206) * temperature - 1930.06;
    ko = (b1 * sr + a1) * salinity + kw;

    dk = (b * pressure + a) * pressure + ko;
    k35 = (5.03217e-5 * pressure + 3.359406) * pressure + 21582.27;
    gam = pressure / k35;
    pk = 1.0 - gam;
    sva = sva * pk + (v350p + sva) * pressure * dk / (k35 * (k35 + dk));
    v350p = v350p * pk;

    dr35p = gam / v350p;
    dvan = sva / (v350p * (v350p + sva));
    *sigma = dr350 + dr35p - dvan;

    return sva * 1.0e+8;
}

double depth(double pressure, double latitude)
{
    double x, gr, depth;

    x = sin(latitude / 57.29578);
    x = x * x;
    gr = 9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * x) * x) + 1.092e-6 *
         pressure;
    depth = (((-1.82e-15 * pressure + 2.279e-10) * pressure - 2.2512e-5) *
             pressure + 9.72659) * pressure;

    return depth / gr;
}

double freezing_point(double salinity, double pressure)
{
    return (-0.0575 + 1.710523e-3 * sqrt(fabs(salinity)) - 2.154996e-4 *
            salinity) * salinity - 7.53e-4 * pressure;
}

double specific_heat(double salinity, double temperature, double pressure)
{
    double a, b, c, cp0, cp1, cp2, sr;

    pressure = pressure / 10.0;
    sr = sqrt(fabs(salinity));

    a = (-1.38385e-3 * temperature + 0.1072763) * temperature - 7.643575;
    b = (5.148e-5 * temperature - 4.07718e-3) * temperature + 0.1770383;
    c = (((2.093236e-5 * temperature - 2.654387e-3) * temperature +
         0.1412855) * temperature -3.720283) * temperature + 4217.4;
    cp0 = (b * sr + a) * salinity + c;

    a = (((1.7168e-8 * temperature + 2.0357e-6) * temperature - 3.13885e-4) *
         temperature + 1.45747e-2) * temperature - 0.49592;
    b = (((2.2956e-11 * temperature - 4.0027e-9) * temperature + 2.87533e-7) *
         temperature - 1.08645e-5) * temperature + 2.4931e-4;
    c = ((6.136e-13 * temperature - 6.5637e-11) * temperature + 2.6380e-9) *
         temperature - 5.422e-8;
    cp1 = ((c * pressure + b) * pressure + a) * pressure;

    a = (((-2.9179e-10 * temperature + 2.5941e-8) * temperature + 9.802e-7) *
         temperature - 1.28315e-4) * temperature + 4.9247e-3;
    b = (3.122e-8 * temperature -1.517e-6) * temperature - 1.2331e-4;
    a = (a + b * sr) * salinity;
    b = ((1.8448e-11 * temperature - 2.3905e-9) * temperature + 1.17054e-7) *
        temperature - 2.9558e-6;
    b = (b + 9.971e-8 * sr) * salinity;
    c = (3.513e-13 * temperature - 1.7682e-11) * temperature + 5.540e-10;
    c = (c - 1.4300e-12 * temperature * sr) * salinity;
    cp2 = ((c * pressure + b) * pressure + a) * pressure;

    return cp0 + cp1 + cp2;
}

double adiabatic_temperature_gradient(double salinity, double temperature,
                                      double pressure)
{
    salinity = salinity - 35.0;

    return (((-2.1687e-16 * temperature + 1.8676e-14) * temperature -
            4.6206e-13) * pressure + ((2.7759e-12 * temperature - 1.1351e-10) *
            salinity + ((-5.4481e-14 * temperature + 8.733e-12) * temperature -
            6.7795e-10) * temperature + 1.8741e-8)) * pressure + (-4.2393e-8 *
            temperature + 1.8932e-6) * salinity + ((6.6228e-10 * temperature -
            6.836e-8) * temperature + 8.5258e-6) * temperature + 3.5803e-5;
}
