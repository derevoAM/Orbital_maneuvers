#ifndef ORBITAL_MANEUVERS_ORBITAL_MANEUVERS_H
#define ORBITAL_MANEUVERS_ORBITAL_MANEUVERS_H

#include "Orbital_elements_convertion.h"

template<typename T>
T Hohmann_transfer(const COE<T> &initial, const COE<T> &final)
{
    T mu = initial.mu;
    T a_trans = (initial.a + final.a) / 2;
    T v_in = std::sqrt(mu / initial.a);
    T v_fin = std::sqrt(mu / final.a);
    T v_trans1 = std::sqrt(2 * mu / initial.a - mu / a_trans);
    T v_trans2 = std::sqrt(2 * mu / final.a - mu / a_trans);
    T delta_v = abs(v_trans1 - v_in) + abs(v_fin - v_trans2);
    return delta_v;

}

template<typename T>
T Bi_elliptic_transfer(const COE<T> &initial, const COE<T> &final, T r_b)
{
    T mu = initial.mu;
    T a_trans1 = (initial.a + r_b) / 2;
    T a_trans2 = (final.a + r_b) / 2;
    T v_in = std::sqrt(mu / initial.a);
    T v_fin = std::sqrt(mu / final.a);
    T v_trans1a = std::sqrt(2 * mu / initial.a - mu / a_trans1);
    T v_trans1b = std::sqrt(2 * mu / r_b - mu / a_trans1);
    T v_trans2b = std::sqrt(2 * mu / r_b - mu / a_trans2);
    T v_trans2c = std::sqrt(2 * mu / final.a - mu / a_trans2);


    T delta_v = abs(v_trans1a - v_in) + abs(v_trans2b - v_trans1b) + abs(v_fin - v_trans2c);
    return delta_v;

}

#endif //ORBITAL_MANEUVERS_ORBITAL_MANEUVERS_H
