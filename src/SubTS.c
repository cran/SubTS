#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <stdlib.h>

void rF1(int *n, double *a, double *p, double *result)
{

    double u1, u2, y;

    for (int i = 0; i < *n; i++)
    {
        for (;;)
        {
            u1 = unif_rand();
            u2 = unif_rand();
            y = pow(1 - log(u2), 1 / *p);
            if (u1 < pow(y, -*p - *a))
            {
                result[i] = y;
                break;
            }
        }
    }
}

void rF2(int *n, double *a, double *p, double *result)
{

    double u1, u2, y;

    for (int i = 0; i < *n; i++)
    {
        for (;;)
        {
            u1 = unif_rand();
            u2 = rbeta((1 - *a) / (*p - 1), 2);
            y = pow(u2, 1 / (*p - 1));
            if (u1 < (exp(-pow(y, *p)) - exp(-y)) / (y - pow(y, *p)))
            {
                result[i] = y;
                break;
            }
        }
    }
}

void simTandW(double *alpha, double *cMax, double *lambda, double *result)
{

    double zeta = 1 / gammafn(1 - *alpha);
    double u, u1, y, au, r, v, u2;

    for (;;)
    {
        u = unif_rand() * M_PI;
        u1 = unif_rand();
        y = 1 - pow(u1, 1 / (1 - *alpha));
        au = pow(sin(*alpha * u), *alpha / (1 - *alpha)) * sin((1 - *alpha) * u) / pow(sin(u), 1 / (1 - *alpha));

        r = rgamma(2 - *alpha, 1 / (au - *lambda));
        v = unif_rand();
        if (v <= au * exp(zeta * pow(r, 1 - *alpha) * pow(y, *alpha) - *lambda * r) * pow(au - *lambda, *alpha - 2) * pow(y, *alpha - 1) * (1 - pow(1 - y, *alpha)) / *cMax)
        {
            u2 = unif_rand();
            result[0] = pow(r, 1 - *alpha) * pow(y, *alpha);
            result[1] = (1 - y) * (pow(1 - u2 * (1 - pow(1 - y, *alpha)), -1 / *alpha) - 1);
            break;
        }
    }
}

void simCondS(double *t, double *alpha, double *result)
{

    double u1, au1, u2, z;

    for (;;)
    {
        u1 = unif_rand() * M_PI;
        au1 = pow(sin(*alpha * u1), *alpha / (1 - *alpha)) * sin((1 - *alpha) * u1) / pow(sin(u1), 1 / (1 - *alpha));
        u2 = unif_rand();
        z = pow(-log(u2) / (au1 * pow(*t, 1 / (1 - *alpha))), 1 - 1 / *alpha);
        if (z < 1)
            break;
    }
    *result = z;
}

void rTrunS(int *n, double *t, double *alpha, double *cMax, double *lambda, double *result)
{

    double z, s, newT;
    double tAndW[] = {0, 0};

    for (int i = 0; i < *n; i++)
    {
        z = 0;
        s = 0;
        for (;;)
        {
            simTandW(alpha, cMax, lambda, tAndW);
            s = s + tAndW[0];
            if (s > *t)
            {
                s = s - tAndW[0];
                break;
            }
            z = z + 1 + tAndW[1];
        }
        newT = *t - s;
        simCondS(&newT, alpha, &(result[i]));
        result[i] += z;
    }
}

void rTrunSCumulative(int *n, double *t, double *alpha, double *cMax, double *lambda, double *result)
{

    double z, s, newT;
    double tAndW[] = {0, 0};
    double *tempResult = (double *)malloc((*n) * sizeof(double));
    if (tempResult == NULL)
    {
        return;
    }

    for (int i = 0; i < *n; i++)
    {
        z = 0;
        s = 0;
        for (;;)
        {
            simTandW(alpha, cMax, lambda, tAndW);
            s = s + tAndW[0];
            if (s > *t)
            {
                s = s - tAndW[0];
                break;
            }
            z = z + 1 + tAndW[1];
        }
        newT = *t - s;
        simCondS(&newT, alpha, &(tempResult[i]));
        result[i] += z + tempResult[i];
    }
    free(tempResult);
}

void rTrunSOptim(int *n, double *t, double *alpha, double *cMax, double *lambda, double *step, double *result)
{

    int tInt = (int)(*t / *step);
    double tFrac = *t - tInt * (*step);

    if (tInt >= 1)
    {
        for (int i = 0; i < tInt; i++)
        {
            rTrunSCumulative(n, step, alpha, cMax, lambda, result);
        }
    }

    if (tFrac > 0)
    {
        rTrunSCumulative(n, &tFrac, alpha, cMax, lambda, result);
    }
}

void rTrunTSCumulative(int *n, double *t, double *alpha, double *cMax, double *lambda, double *mu, double *result)
{

    double v;
    int temp = 1;
    double *tempResult = (double *)malloc((*n) * sizeof(double));
    if (tempResult == NULL)
    {
        return;
    }
    for (int i = 0; i < *n; i++)
    {
        for (;;)
        {
            rTrunS(&temp, t, alpha, cMax, lambda, &(tempResult[i]));
            v = unif_rand();
            if (v <= exp(-*mu * tempResult[i]))
            {
                result[i] += tempResult[i];
                break;
            }
        }
    }
    free(tempResult);
}

void rTrunTS(int *n, double *t, double *alpha, double *cMax, double *lambda, double *mu, double *step, double *result)
{

    int tInt = (int)(*t / *step);
    double tFrac = *t - tInt * (*step);

    if (tInt >= 1)
    {
        for (int i = 0; i < tInt; i++)
        {
            rTrunTSCumulative(n, step, alpha, cMax, lambda, mu, result);
        }
    }

    if (tFrac > 0)
    {
        rTrunTSCumulative(n, &tFrac, alpha, cMax, lambda, mu, result);
    }
}

void rPRDTS(int *n, double *t, double *alpha, double *cMax, double *lambda, double *p, double *k1, double *k2, double *step, double *result)
{

    int n1, n2;
    double cp1, cp2, temp;
    int one = 1;
    double oneD = 1;
    rTrunTS(n, t, alpha, cMax, lambda, &oneD, step, result);

    for (int i = 0; i < *n; i++)
    {
        cp1 = 0;
        cp2 = 0;
        n1 = rpois(*t * *k1 * *alpha / gammafn(1 - *alpha));
        n2 = rpois(*t * *k2 * *alpha / gammafn(1 - *alpha));

        for (int j = 0; j < n1; j++)
        {
            rF1(&one, alpha, p, &temp);
            cp1 += temp;
        }

        for (int j = 0; j < n2; j++)
        {
            rF2(&one, alpha, p, &temp);
            cp2 += temp;
        }
        result[i] += (cp1 + cp2);
    }
}

void rPRDTSneg(int *n, double *t, double *alpha, double *p, double *result)
{

    double k = gammafn(-*alpha / *p) / *p;
    int n1;
    double cp;

    for (int i = 0; i < *n; i++)
    {
        cp = 0;
        n1 = rpois(*t * k);

        for (int j = 0; j < n1; j++)
        {
            cp += pow(rgamma(-*alpha / *p, 1), 1 / *p);
        }
        result[i] = cp;
    }
}

void rDickman(int *n, double *t, double *result)
{

    double s, g, y, v, ct, m, tt;
    g = -digamma(1);

    for (int i = 0; i < *n; i++)
    {
        s = 0;
        tt = *t;
        for (;;)
        {
            for (;;)
            {
                ct = -1.25 * log(unif_rand());
                y = rbeta(ct, 0.5);
                m = y - 1 + pow(1 - y, 1 - unif_rand());
                v = unif_rand();
                if (v <= (25.0 / 47.0) * gammafn(0.5) * exp((0.8 - g) * ct) * (-log(1 - y)) * sqrt(1 - y) / gammafn(ct + 0.5))
                {
                    break;
                }
            }
            if (ct > tt)
            {
                break;
            }
            s += (1 + m);
            tt -= ct;
        }
        result[i] = s + pow(unif_rand(), 1 / tt);
    }
}

void rTrunGammaCumulative(int *n, double *t, double *mu, double *result)
{

    double v;
    int temp = 1;
    double *tempResult = (double *)malloc((*n) * sizeof(double));
    if (tempResult == NULL)
    {
        return;
    }
    for (int i = 0; i < *n; i++)
    {
        for (;;)
        {
            rDickman(&temp, t, &(tempResult[i]));
            v = unif_rand();
            if (v <= exp(-*mu * tempResult[i]))
            {
                result[i] += tempResult[i];
                break;
            }
        }
    }
    free(tempResult);
}

void rTrunGamma(int *n, double *t, double *mu, double *step, double *result)
{

    int tInt = (int)(*t / *step);
    double tFrac = *t - tInt * (*step);

    if (tInt >= 1)
    {
        for (int i = 0; i < tInt; i++)
        {
            rTrunGammaCumulative(n, step, mu, result);
        }
    }

    if (tFrac > 0)
    {
        rTrunGammaCumulative(n, &tFrac, mu, result);
    }
}

void rPGamma(int *n, double *t, double *p, double *k1, double *k2, double *step, double *result)
{

    int n1, n2;
    double cp1, cp2, temp;
    int one = 1;
    double oneD = 1;
    double zero = 0;
    rTrunGamma(n, t, &oneD, step, result);

    for (int i = 0; i < *n; i++)
    {
        cp1 = 0;
        cp2 = 0;
        n1 = rpois(*t * *k1);
        n2 = rpois(*t * *k2);

        for (int j = 0; j < n1; j++)
        {
            rF1(&one, &zero, p, &temp);
            cp1 += temp;
        }

        for (int j = 0; j < n2; j++)
        {
            rF2(&one, &zero, p, &temp);
            cp2 += temp;
        }
        result[i] += (cp1 + cp2);
    }
}

double vFunction(double alpha, double theta, double theta0)
{
    if (alpha == 1)
    {
        return (1 - 2 * theta / M_PI) * exp((theta - M_PI / 2) * tan(theta)) / cos(theta);
    }
    else
    {
        return pow(cos(alpha * theta0) * cos(theta), 1 / (alpha - 1)) * cos((alpha - 1) * theta + alpha * theta0) / pow(fabs(sin(alpha * (theta + theta0))), alpha / (alpha - 1));
    }
}

double muFunction(double alpha, double ell, double epsilon, double theta, double theta0)
{
    double v = vFunction(alpha, theta, theta0);
    if (alpha == 1)
    {
        return 2 * log(2 / (epsilon * ell * M_PI * v)) / M_PI;
    }
    else
    {
        return pow((alpha - 1) / (ell * epsilon * alpha * v), alpha - 1);
    }
}

double zetaFunction(double alpha, double ell, double epsilon, double theta, double theta0)
{
    if (alpha == 1)
    {
        return M_PI / (2 * ell);
    }
    else
    {
        double v = vFunction(alpha, theta, theta0);

        return pow(v * epsilon * alpha * ell, alpha - 1) / (ell * pow(alpha - 1, alpha));
    }
}

void simH1(double alpha, double theta0, double *result)
{
    double theta = (unif_rand() - 0.5) * M_PI;
    double w = -log(unif_rand());
    double x;
    if (alpha == 1)
    {
        x = 2 / M_PI * ((M_PI / 2 - theta) * tan(theta) + log((M_PI * w * cos(theta)) / (M_PI - 2 * theta)));
    }
    else
    {
        x = pow(vFunction(alpha, theta, theta0) / w, 1 / alpha - 1);
        if (theta0 + theta < 0)
            x = (-1) * x;
    }

    result[0] = x;
    result[1] = theta;
}

void simH2(double alpha, double ell, double epsilon, double theta0, double *result)
{
    double z = norm_rand();
    double theta = unif_rand() * (M_PI / 2 + theta0) - theta0;
    double x = muFunction(alpha, ell, epsilon, theta, theta0) + fabs(z) / sqrt(zetaFunction(alpha, ell, epsilon, theta, theta0));
    result[0] = x;
    result[1] = theta;
}

void simH(double alpha, double ell, double epsilon, double p, double theta0, double *result)
{
    double u = unif_rand();
    if (u <= p)
    {
        simH1(alpha, theta0, result);
    }
    else
    {
        simH2(alpha, ell, epsilon, theta0, result);
    }
}

double simGstar(double alpha, double ell, double epsilon, double p, double theta0)
{
    double temp[] = {0, 0};
    double u, x, theta, phi, k, g, h, h1, h2, v, mu, zeta;
    while (true)
    {
        u = unif_rand();
        simH(alpha, ell, epsilon, p, theta0, temp);
        x = temp[0];
        theta = temp[1];
        v = vFunction(alpha, theta, theta0);
        mu = muFunction(alpha, ell, epsilon, theta, theta0);
        zeta = zetaFunction(alpha, ell, epsilon, theta, theta0);

        if ((x <= mu) || (theta <= -theta0))
        {
            h2 = 0;
        }
        else
        {
            h2 = alpha * exp(-zeta * pow(x - mu, 2) / 2) * sqrt(2 * zeta / M_PI) / M_PI;
        }

        if (alpha == 1)
        {
            // This is K
            k = fmax(exp(2 * (1 - log(epsilon)) / (M_PI * ell)) / p, M_PI * sqrt(ell) * exp(-2 / (M_PI * ell) - 1) / (2 * (1 - epsilon) * (1 - p)));

            // This is f^*(x,theta)
            h1 = .5 * exp(x * M_PI / 2 - exp(M_PI * x / 2) * v) * v;

            // This is g^*
            g = pow(ell, 2 / (M_PI * ell)) * exp(x / ell) * h1;
        }
        else
        {
            // this is K/C
            k = fmax(exp(fmax(1, pow(ell * epsilon, 1 - alpha) * alpha / fabs(cos(M_PI * alpha / 2))) / ell) / p, exp(-pow(epsilon, 1 - alpha) / (pow(ell, alpha) * cos(M_PI * alpha / 2)) + (alpha - 3) / 2) * sqrt(M_PI) * pow(epsilon * alpha, (1 - alpha) / 2) * pow(ell / (alpha - 1), 1 - alpha / 2) * pow((3 - alpha) / (1 - epsilon), (3 - alpha) / 2) / pow(2, 2 - alpha / 2) / (1 - p));

            // This is f^*(x,theta)
            if ((x < 0 && theta >= -theta0) || (x > 0 && theta <= -theta0))
            {
                h1 = 0;
            }
            else
            {
                h1 = exp(-pow(fabs(x), alpha / (alpha - 1)) * v) * alpha * v * pow(fabs(x), 1 / (alpha - 1)) / (M_PI * (alpha - 1));
            }

            // this is g^*/C
            g = exp(x / ell) * h1;
        }
        h = p * h1 + (1 - p) * h2;
        phi = g / (k * h);
        if (u <= phi)
            return x;
    }
}

void rCTSM(int *n, double *alpha, double *c, double *ell, double *epsilon, double *p, double *result)
{
    double shift;
    if (*alpha == 1)
    {
        shift = 2 * (*c) * log(*c) / M_PI;
    }
    else
    {
        shift = 0;
    }

    double theta0 = M_PI * (1 / (*alpha) - 0.5);
    for (int i = 0; i < *n; i++)
    {
        result[i] = simGstar(*alpha, *ell, *epsilon, *p, theta0) * (*c) - shift;
    }
}
