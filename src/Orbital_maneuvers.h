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

#endif //ORBITAL_MANEUVERS_ORBITAL_MANEUVERS_H
