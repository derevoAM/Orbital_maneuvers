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
    //std::cout << r2 << "\n";
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

    T delta_v = v1h + v2h - std::sqrt(v0r * v0r + (v0t + v1ha) * (v0t + v1ha)) -
                std::sqrt(v3r * v3r + (v3t - v2ha) * (v3t - v2ha));
    //std::cout << v3r << " " << v3t << "\n";
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
    //std::cout << v3r << " " << v3t << " " << norm(r2) << "\n";

    T delta_v = std::sqrt(v3r * v3r + (v3t + v1ht) * (v3t + v1ht)) - std::sqrt(v0r * v0r + (v0t + v2ht) * (v0t + v2ht));
    return delta_v;
}

template<typename T>
T Inclination_only_transfer(const COE<T> &initial, const COE<T> &final) {
    T r1 = initial.p / (1 + initial.e * cos(2 * M_PI - initial.w));
    T r2 = initial.p / (1 + initial.e * cos(M_PI - initial.w));

    T v1 = std::sqrt(2 * initial.mu / r1 - initial.mu / initial.a);
    T v2 = std::sqrt(2 * initial.mu / r2 - initial.mu / initial.a);

    T fi1 = atan(initial.e * sin(2 * M_PI - initial.w) / (1 + initial.e * cos(2 * M_PI - initial.w)));
    T fi2 = atan(initial.e * sin(M_PI - initial.w) / (1 + initial.e * cos(M_PI - initial.w)));

    T delta_v1 = 2 * v1 * cos(fi1) * sin(abs(initial.i - final.i) / 2);
    T delta_v2 = 2 * v2 * cos(fi2) * sin(abs(initial.i - final.i) / 2);
    return std::min(delta_v1, delta_v2);
}

template<typename T>
std::tuple<T, std::vector<T>, std::vector<T>> General_plane_change(COE<T> &initial, const COE<T> &final) {
    auto [r_i, v_i] = COE2RV(initial);
    auto [r_f, v_f] = COE2RV(final);
    T mu = initial.mu;
    T delta_v1, delta_v2;

    // finding normal vectors
    std::vector<T> h1 = cross_product(r_i, v_i); // coordinates are A, B, C of plane equation
    h1 = h1 / norm(h1);

    std::vector<T> h2 = cross_product(r_f, v_f);
    h2 = h2 / norm(h2);

    std::vector<T> a = cross_product(h1, h2); // vector of plane intersection
    a = a / norm(a);
    std::cout << "a= " << a << "\n";

    std::vector<T> e = (r_i * (scalar(v_i, v_i) - mu / norm(r_i)) - v_i * scalar(r_i, v_i)) / mu;
    //std::cout << e << "\n";


    T nu = acos(scalar(e, a) / (norm(e) * norm(a)));
    if (scalar(a, v_i) < 0) nu = 2 * M_PI - nu;

    std::cout << "nu = " << nu * 180 / M_PI << "\n";

    initial.nu = nu;
    std::vector<T> v2_1, v2_2; //two vectors for 2 nodes

    // 1 node of intersecting planes
    auto [r11, v11] = COE2RV(initial);
    T alpha = acos(scalar(h1, h2));
    delta_v1 = std::sqrt(2 * scalar(v11, v11) * (1 - cos(alpha)));

    std::cout << v11 << "\n";
    if (scalar(h1, h2) >= 1e-5) {
        v2_1 = v11 * cos(alpha) + h1 * norm(v11) * sin(alpha);
    } else {
        v2_1 = v11 * (-cos(alpha)) + h1 * norm(v11) * sin(M_PI - alpha);
        //if (scalar(h2, v2_1) / norm(v2_1) > 1e-1) v2_1 = v11 * (-cos(alpha)) - h1 * norm(v11) * sin(M_PI - alpha);
    }
    std::cout << v2_1 << "\n";
    std::cout << "v2=v1 " << norm(v2_1) / norm(v11) << "\n";

    // 2 node of intersecting planes
    if (nu < M_PI) initial.nu = nu + M_PI;
    else initial.nu = nu - M_PI;

    std::cout << "nu = " << initial.nu * 180 / M_PI << "\n";

    auto [r12, v12] = COE2RV(initial);
    delta_v2 = std::sqrt(2 * scalar(v12, v12) * (1 - cos(alpha)));

    if (scalar(h1, h2) >= 1e-30) {
        v2_2 = v12 * cos(alpha) - h1 * norm(v12) * sin(alpha);
        //std::cout << "haha";
    } else {
        v2_2 = v12 * (-cos(alpha)) + h1 * norm(v12) * sin(M_PI - alpha);
        //if (scalar(h2, v2_2) > 1e-30) v2_2 = v12 * (-cos(alpha)) - h1 * norm(v12) * sin(M_PI - alpha);
    }
    //std::cout << "v2=v1 " << norm(v2_2) / norm(v12) << "\n";
    //std::cout << r11 << "\n" << r12 << "\n";
    //std::cout << v12 << "\n" << v2_2 << "\n";

    std::cout << delta_v1 << " " << delta_v2 << "\n";
//    if (delta_v1 < delta_v2) return {delta_v1, r11, v2_1};
//    else return {delta_v2, r12, v2_2};
    return {delta_v1, r11, v2_1};
}


#endif //ORBITAL_MANEUVERS_ORBITAL_MANEUVERS_H
