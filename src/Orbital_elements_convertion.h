#ifndef ORBITAL_MANEUVERS_ORBITAL_ELEMENTS_CONVERTION_H
#define ORBITAL_MANEUVERS_ORBITAL_ELEMENTS_CONVERTION_H

#include <iostream>
#include <math.h>
#include "Vector.h"
#include "Dense.h"


template<typename T>
struct COE {
    T p;
    T a;
    T e;
    T i;
    T W; //Omega big
    T w; //omega small
    T nu;
    T u;
    T lam_true;
    T w_true;
    T mu;
    int flag = 0; // determine which type of orbit: 1 - cir eq, 2 - cir in, 3 - el eq, 4 - el in
};

template<typename T>
COE<T> RV2COE(const std::vector<T> &r, const std::vector<T> &v, T mu) {
    std::vector<T> h = cross_product(r, v);
    std::vector<T> n = cross_product(std::vector<T>{0, 0, 1}, h);
    std::vector<T> e = (r * (scalar(v, v) - mu / norm(r)) - v * scalar(r, v)) / mu;
    T ksi = scalar(v, v) / 2 - mu / norm(r);
    T a = -mu / (2 * ksi);
    T p = scalar(h, h) / mu;
    T i = acos(h[2] / norm(h));
    T w_true, W, w, lam_true, u, nu;
    int flag = 0;

    if (norm(e) < 1e-4) { // circular case
        nu = acos(
                r[0] / norm(r)); // let's consider, that in circular case true anomaly is the angle between I axis and r
        if (scalar(r, v) < 0) nu = 2 * M_PI - nu;
        if (norm(n) == 0) // circular equatorial case
        {
            lam_true = r[0] / norm(r);
            if (r[1] < 0) lam_true = 2 * M_PI - lam_true;

            // if degree parameter is equal 10, then it's not defined
            W = 10;
            w = 10;
            u = 10;
            w_true = 10;
            flag = 1;

        } else // circular inclined case
        {
            u = acos(scalar(n, r) / (norm(n) * norm(r)));
            if (r[2] < 0) u = 2 * M_PI - u;
            W = acos(n[0] / norm(n));
            if (n[1] < 0) W = 2 * M_PI - W;
            w = 10;
            w_true = 10;
            lam_true = 10;
            flag = 2;
        }
    } else // elliptic case
    {
        nu = acos(scalar(e, r) / (norm(e) * norm(r)));
        if (scalar(r, v) < 0) nu = 2 * M_PI - nu;
        if (norm(n) == 0) // elliptic equatorial case
        {
            w_true = e[0] / norm(e);
            if (e[1] < 0) w_true = 2 * M_PI - w_true;
            W = 10;
            w = 10;
            lam_true = 10;
            u = 10;
            flag = 3;


        } else // elliptic inclined case
        {
            W = acos(n[0] / norm(n));
            if (n[1] < 0) W = 2 * M_PI - W;
            w = acos(scalar(n, e) / (norm(n) * norm(e)));
            if (e[2] < 0) w = 2 * M_PI - w;
            u = 10;
            w_true = 10;
            lam_true = 10;
            flag = 4;
        }
    }
    return COE<T>{p, a, norm(e), i, W, w, nu, u, lam_true, w_true, mu, flag};
}

template<typename T>
std::pair<std::vector<T>, std::vector<T>> COE2RV(const COE<T> &elem) {
    T w, W, nu, i = elem.i;
    if (elem.flag == 1) {
        w = 0;
        W = 0;
        nu = elem.lam_true;
    }
    if (elem.flag == 2) {
        w = 0;
        W = elem.W;
        nu = elem.u;
    }
    if (elem.flag == 3) {
        W = 0;
        w = elem.w_true;
        nu = elem.nu;
    }
    if (elem.flag == 4) {
        nu = elem.nu;
        w = elem.w;
        W = elem.W;
    }
    std::vector<T> r_pqw{elem.p * cos(nu) / (1 + elem.e * cos(nu)), elem.p * sin(nu) / (1 + elem.e * cos(nu)), 0};
    std::vector<T> v_pqw{-std::sqrt(elem.mu / elem.p) * sin(nu), std::sqrt(elem.mu / elem.p) * (elem.e + cos(nu)), 0};
    Dense<T> A(3, 3, {cos(W) * cos(w) - sin(W) * sin(w) * cos(i), -cos(W) * sin(w) - sin(W) * cos(w) * cos(i), sin(W) * sin(i),
                      sin(W) * cos(w) + cos(W) * sin(w) * cos(i), -sin(W) * sin(w) + cos(W) * cos(w) * cos(i), -cos(W) * sin(i),
                      sin(w) * sin(i),                             cos(w) * sin(i),                             cos(i)});
    return std::pair(A * r_pqw, A * v_pqw);
}

#endif //ORBITAL_MANEUVERS_ORBITAL_ELEMENTS_CONVERTION_H
