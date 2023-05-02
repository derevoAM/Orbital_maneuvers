#ifndef ORBITAL_MANEUVERS_ORBITAL_MANEUVERS_H
#define ORBITAL_MANEUVERS_ORBITAL_MANEUVERS_H

#include "Orbital_elements_convertion.h"

template<typename T>
T Hohmann_transfer(const COE<T> &initial, const COE<T> &final) {
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
T Bi_elliptic_transfer_circular_orbits(const COE<T> &initial, const COE<T> &final, T r_b) {
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

template<typename T>
T Bi_elliptic_transfer_elliptic_orbits(const COE<T> &initial, const COE<T> &final, T r_a) {
    auto [r1, v1] = COE2RV(initial);
    auto [r2, v2] = COE2RV(final);
    T mu = initial.mu;

    T p1h = 2 * norm(r1) * r_a / (norm(r1) + r_a);
    T v1h = std::sqrt(mu * p1h) / norm(r1);

    T p2h = 2 * norm(r2) * r_a / (norm(r2) + r_a);
    T v2h = std::sqrt(mu * p2h) / norm(r2);

    T v0r = scalar(v1, r1) / norm(r1);
    T v0t = std::sqrt(scalar(v1, v1) - v0r * v0r);

    T v3r = scalar(v2, r2) / norm(r2);
    T v3t = std::sqrt(scalar(v2, v2) - v3r * v3r);


    T v1ha = std::sqrt(mu * p1h) / r_a;
    T v2ha = std::sqrt(mu * p2h) / r_a;

    std::cout << v3r << " " << v3t << "\n";

    T delta_v = v1h + v2h - std::sqrt(v0r * v0r + (v0t + v1ha) * (v0t + v1ha)) -
                std::sqrt(v3r * v3r + (v3t - v2ha) * (v3t - v2ha));
    return delta_v;

}

template<typename T>
T Two_impulse_transfer_elliptic_orbits(const COE<T> &initial, const COE<T> &final) {
    auto [r1, v1] = COE2RV(initial);
    auto [r2, v2] = COE2RV(final);
    T mu = initial.mu;

    T p_h = 2 * norm(r1) * norm(r2) / (norm(r1) + norm(r2));
    T v1ht = std::sqrt(mu * p_h) / norm(r1);
    T v2ht = std::sqrt(mu * p_h) / norm(r2);

    T v0r = scalar(v1, r1) / norm(r1);
    T v0t = std::sqrt(scalar(v1, v1) - v0r * v0r);

    T v3r = scalar(v2, r2) / norm(r2);
    T v3t = std::sqrt(scalar(v2, v2) - v3r * v3r);
    std::cout << v3r << " " << v3t << " " << norm(r2) << "\n";

    T delta_v = std::sqrt(v3r * v3r + (v3t + v1ht) * (v3t + v1ht)) - std::sqrt(v0r * v0r + (v0t + v2ht) * (v0t + v2ht));
    return delta_v;
}


#endif //ORBITAL_MANEUVERS_ORBITAL_MANEUVERS_H
