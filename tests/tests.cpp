#include "gtest/gtest.h"
#include "../src/Orbital_elements_convertion.h"
#include "../src/Orbital_maneuvers.h"

TEST(ORBITAL_MANEUVERS, RV2COE_ELLIPTIC_INCLINED) {
    /// Converting RV vectors to Keplerian elements for elliptic inclined orbit ///

    std::vector<double> r{6524.834, 6862.875, 6448.296};
    std::vector<double> v{4.901327, 5.533756, -1.976341};
    double mu = 398600.4415;
    COE<double> elem = RV2COE(r, v, mu);
    ASSERT_NEAR(elem.p, 11067.790, 1e-2);
    ASSERT_NEAR(elem.a, 36127.343, 1e-2);
    ASSERT_NEAR(elem.e, 0.832853, 1e-2);
    ASSERT_NEAR(elem.i, 87.87 * M_PI / 180, 1e-2);
    ASSERT_NEAR(elem.W, 227.898 * M_PI / 180, 1e-2);
    ASSERT_NEAR(elem.w, 53.38 * M_PI / 180, 1e-2);
    ASSERT_EQ(elem.flag, 4);

}

TEST(ORBITAL_MANEUVERS, RV2COE_ELLIPTIC_EQUATORIAL) {
    /// Converting RV vectors to Keplerian elements for elliptic equatorial orbit ///

    std::vector<double> r{6524.834, 6862.875, 0};
    std::vector<double> v{4.901327, 1.533756, 0};
    double mu = 398600.4415;
    COE<double> elem = RV2COE(r, v, mu);
    ASSERT_NEAR(elem.a, 6894.965, 1e-1);
    ASSERT_NEAR(elem.e, 0.8926, 1e-2);
    ASSERT_NEAR(elem.i, 180 * M_PI / 180, 1e-2);
    ASSERT_NEAR(elem.w_true, 197.848 * M_PI / 180, 2e-1);
    ASSERT_EQ(elem.flag, 3);

}

TEST(ORBITAL_MANEUVERS, RV2COE_CIRCULAR_INCLINED) {
    /// Converting RV vectors to Keplerian elements for circular inclined orbit ///

    std::vector<double> r{5011.173, 6234.5, -138.135};
    std::vector<double> v{5.499, -4.422, -0.1048};
    double mu = 398600.4415;
    COE<double> elem = RV2COE(r, v, mu);
    ASSERT_NEAR(elem.a, 8000, 5);
    ASSERT_NEAR(elem.e, 0, 1e-2);
    ASSERT_NEAR(elem.i, 178.7 * M_PI / 180, 1e-1);
    ASSERT_NEAR(elem.W, 280.5 * M_PI / 180, 1e-2);
    ASSERT_EQ(elem.flag, 2);

}

TEST(ORBITAL_MANEUVERS, RV2COE_CIRCULAR_EQUATORIAL) {
    /// Converting RV vectors to Keplerian elements for circular equatorial orbit ///

    std::vector<double> r{5011.173, 6234.5, 0};
    std::vector<double> v{5.499, -4.422, 0};
    double mu = 398600.4415;
    COE<double> elem = RV2COE(r, v, mu);
    ASSERT_NEAR(elem.a, 7997, 5);
    ASSERT_NEAR(elem.e, 0, 1e-2);
    ASSERT_NEAR(elem.i, 180 * M_PI / 180, 1e-2);
    ASSERT_EQ(elem.flag, 1);

}



TEST(ORBITAL_MANEUVERS, COE2RV_ELLIPTIC_INCLINED) {
    COE<double> elem;
    elem.p = 11067.790;
    elem.e = 0.83285;
    elem.i = 87.87 * M_PI / 180;
    elem.W = 227.898 * M_PI / 180;
    elem.w = 53.38 * M_PI / 180;
    elem.nu = 92.335 * M_PI / 180;
    elem.flag = 4;
    elem.mu = 398600.4415;
    auto [r, v] = COE2RV(elem);
    ASSERT_NEAR(r[0], 6525.344, 1);
    ASSERT_NEAR(r[1], 6861.535, 1);
    ASSERT_NEAR(r[2], 6449.125, 1);
    ASSERT_NEAR(v[0], 4.902276, 1e-3);
    ASSERT_NEAR(v[1], 5.533124, 1e-3);
    ASSERT_NEAR(v[2], -1.975709, 1e-3);

}

TEST(ORBITAL_MANEUVERS, HOHMANN_TRANSFER) {
    COE<double> elem1;
    elem1.a = 191.3441 + 6378.137;
    COE<double> elem2;
    elem2.a = 35781.34857 + 6378.137;
    elem1.mu = 398600.4415;
    elem2.mu = 398600.4415;
    double delta_v = Hohmann_transfer(elem1, elem2);
    ASSERT_NEAR(delta_v, 3.935224, 1e-6);

}

TEST(ORBITAL_MANEUVERS, BI_ELLIPTIC_TRANSFER_CIRCULAR_ORBITS) {
    COE<double> elem1;
    elem1.a = 191.3441 + 6378.137;
    COE<double> elem2;
    elem2.a = 376310 + 6378.137;
    elem1.mu = 398600.4415;
    elem2.mu = 398600.4415;
    double r_b = 503873 + 6378.137;
    double delta_v = Bi_elliptic_transfer_circular_orbits(elem1, elem2, r_b);
    ASSERT_NEAR(delta_v, 3.904057, 1e-6);

}

TEST(ORBITAL_MANEUVERS, BI_ELLIPTIC_TRANSFER_BOUNDARY_CONDITION_EXAMPLE_1) {
    COE<double> elem1;
    elem1.nu = 10 * M_PI / 180;
    elem1.e = 0.2;
    elem1.p = 10320 * (1 - elem1.e);
    elem1.a = 10320 / (1 + elem1.e);
    elem1.i = 0;
    elem1.flag = 3;
    elem1.w_true = M_PI / 2;
    elem1.mu = 398600.4415;

    COE<double> elem2;
    elem2.nu = 0 * M_PI / 180;
    elem2.e = 0.2;
    elem2.p = 10320 * 13.43 * (1 - elem2.e);
    elem2.a = 10320 * 13.43 / (1 + elem2.e);
    elem2.i = 0;
    elem2.flag = 3;
    elem2.w_true = 35 * M_PI / 180;
    elem2.mu = 398600.4415;
    double delta_v = Bi_elliptic_transfer_elliptic_orbits(elem1, elem2, 10320 * 13.43);
    ASSERT_NEAR(delta_v, 3.16723, 1e-5);
}

TEST(ORBITAL_MANEUVERS, BI_ELLIPTIC_TRANSFER_GENERAL_EXAMPLE_1) {
    COE<double> elem1;
    elem1.nu = 10 * M_PI / 180;
    elem1.e = 0.2;
    elem1.p = 10320 * (1 - elem1.e);
    elem1.a = 10320 / (1 + elem1.e);
    elem1.i = 0;
    elem1.flag = 3;
    elem1.w_true = M_PI / 2;
    elem1.mu = 398600.4415;

    COE<double> elem2;
    elem2.nu = 0 * M_PI / 180;
    elem2.e = 0.2;
    elem2.p = 10320 * 13.43 * (1 - elem2.e);
    elem2.a = 10320 * 13.43 / (1 + elem2.e);
    elem2.i = 0;
    elem2.flag = 3;
    elem2.w_true = 35 * M_PI / 180;
    elem2.mu = 398600.4415;
    double delta_v = Bi_elliptic_transfer_elliptic_orbits(elem1, elem2, 2 * 10320 * 13.43);
    ASSERT_NEAR(delta_v, 3.15139, 1e-4);
}

TEST(ORBITAL_MANEUVERS, TWO_IMPULSE_TRANSFER_EXAMPLE_1) {
    COE<double> elem1;
    elem1.nu = 10 * M_PI / 180;
    elem1.e = 0.2;
    elem1.p = 10320 * (1 - elem1.e);
    elem1.a = 10320 / (1 + elem1.e);
    elem1.i = 0;
    elem1.flag = 3;
    elem1.w_true = M_PI / 2;
    elem1.mu = 398600.4415;

    COE<double> elem2;
    elem2.nu = 0 * M_PI / 180;
    elem2.e = 0.2;
    elem2.p = 10320 * 13.43 * (1 - elem2.e);
    elem2.a = 10320 * 13.43 / (1 + elem2.e);
    elem2.i = 0;
    elem2.flag = 3;
    elem2.w_true = 35 * M_PI / 180;
    elem2.mu = 398600.4415;
    double delta_v = Two_impulse_transfer_elliptic_orbits(elem1, elem2);
    ASSERT_NEAR(delta_v, 3.55158, 1e-4);
}

TEST(ORBITAL_MANEUVERS, BI_ELLIPTIC_TRANSFER_BOUNDARY_CONDITION_EXAMPLE_2) {
    COE<double> elem1;
    elem1.nu = 10 * M_PI / 180;
    elem1.e = 0.2;
    elem1.p = 10320 * (1 - elem1.e);
    elem1.a = 10320 / (1 + elem1.e);
    elem1.i = 30 * M_PI / 180;
    elem1.flag = 4;
    elem1.w = 15 * M_PI / 180;
    elem1.W = 45 * M_PI / 180;
    elem1.mu = 398600.4415;

    COE<double> elem2;
    elem2.nu = 50 * M_PI / 180;
    elem2.e = 0.2;
    elem2.p = 10320 * 18.98 * (1 - elem1.e);
    elem2.a = 10320 * 18.98 / (1 + elem1.e);
    elem2.i = 30 * M_PI / 180;
    elem2.flag = 4;
    elem2.w = 15 * M_PI / 180;
    elem2.W = 45 * M_PI / 180;
    elem2.mu = 398600.4415;
    double delta_v = Bi_elliptic_transfer_elliptic_orbits(elem1, elem2, 10320 * 18.98);
    ASSERT_NEAR(delta_v, 3.14986, 1e-4);
}

TEST(ORBITAL_MANEUVERS, BI_ELLIPTIC_TRANSFER_GENERAL_EXAMPLE_2) {
    COE<double> elem1;
    elem1.nu = 10 * M_PI / 180;
    elem1.e = 0.2;
    elem1.p = 10320 * (1 - elem1.e);
    elem1.a = 10320 / (1 + elem1.e);
    elem1.i = 30 * M_PI / 180;
    elem1.flag = 4;
    elem1.w = 15 * M_PI / 180;
    elem1.W = 45 * M_PI / 180;
    elem1.mu = 398600.4415;

    COE<double> elem2;
    elem2.nu = 50 * M_PI / 180;
    elem2.e = 0.2;
    elem2.p = 10320 * 18.98 * (1 - elem1.e);
    elem2.a = 10320 * 18.98 / (1 + elem1.e);
    elem2.i = 30 * M_PI / 180;
    elem2.flag = 4;
    elem2.w = 15 * M_PI / 180;
    elem2.W = 45 * M_PI / 180;
    elem2.mu = 398600.4415;
    double delta_v = Bi_elliptic_transfer_elliptic_orbits(elem1, elem2, 2 * 10320 * 18.98);
    ASSERT_NEAR(delta_v, 3.11046, 1e-4);
}

TEST(ORBITAL_MANEUVERS, TWO_IMPULSE_TRANSFER_EXAMPLE_2) {
    COE<double> elem1;
    elem1.nu = 10 * M_PI / 180;
    elem1.e = 0.2;
    elem1.p = 10320 * (1 - elem1.e);
    elem1.a = 10320 / (1 + elem1.e);
    elem1.i = 30 * M_PI / 180;
    elem1.flag = 4;
    elem1.w = 15 * M_PI / 180;
    elem1.W = 45 * M_PI / 180;
    elem1.mu = 398600.4415;

    COE<double> elem2;
    elem2.nu = 50 * M_PI / 180;
    elem2.e = 0.2;
    elem2.p = 10320 * 18.98 * (1 - elem1.e);
    elem2.a = 10320 * 18.98 / (1 + elem1.e);
    elem2.i = 30 * M_PI / 180;
    elem2.flag = 4;
    elem2.w = 15 * M_PI / 180;
    elem2.W = 45 * M_PI / 180;
    elem2.mu = 398600.4415;
    double delta_v = Two_impulse_transfer_elliptic_orbits(elem1, elem2);
    ASSERT_NEAR(delta_v, 3.45414, 1e-4);
}

TEST(ORBITAL_MANEUVERS, BI_ELLIPTIC_TRANSFER_CIRCULAR_ORBITS_GENERAL) {
    COE<double> elem1;
    elem1.a = 191.3441 + 6378.137;
    elem1.p = elem1.a;
    COE<double> elem2;
    elem2.a = 376310 + 6378.137;
    elem2.p = elem2.a;
    elem1.mu = 398600.4415;
    elem2.mu = 398600.4415;

    elem1.flag = 1;
    elem2.flag = 1;
    elem1.i = 0;
    elem2.i = 0;
    elem1.lam_true = 10 * M_PI / 180;
    elem2.lam_true = 10 * M_PI / 180;
    double r_b = 503873 + 6378.137;
    double delta_v = Bi_elliptic_transfer_elliptic_orbits(elem1, elem2, r_b);
    ASSERT_NEAR(delta_v, 3.904057, 1e-6);

}

TEST(ORBITAL_MANEUVERS, INCLINATION_ONLY_TRANSFER) {
    COE<double> elem1;
    elem1.e = 0.3;
    elem1.p = 17858.7836;
    elem1.a = elem1.p / (1 - elem1.e * elem1.e);
    elem1.i = 30 * M_PI / 180;
    elem1.flag = 4;
    elem1.w = 30 * M_PI / 180;
    elem1.W = 45 * M_PI / 180;
    elem1.mu = 398600.4415;

    COE<double> elem2;
    elem2.e = 0.3;
    elem2.p = 17858.7836;
    elem2.a = elem1.p / (1 - elem1.e * elem1.e);
    elem2.i = 45 * M_PI / 180;
    elem2.flag = 4;
    elem2.w = 30 * M_PI / 180;
    elem2.W = 45 * M_PI / 180;
    elem2.mu = 398600.4415;
    double delta_v = Inclination_only_transfer(elem1, elem2);
    ASSERT_NEAR(delta_v, 0.912833, 1e-4);
}

TEST(ORBITAL_MANEUVERS, PLANE_CHANGE) {
    COE<double> elem1;
    elem1.e = 0.3;
    elem1.p = 17858.7836;
    elem1.a = elem1.p / (1 - elem1.e * elem1.e);
    elem1.i = 70 * M_PI / 180;
    elem1.flag = 4;
    elem1.w = 30 * M_PI / 180;
    elem1.W = 90 * M_PI / 180;
    elem1.mu = 398600.4415;

    COE<double> elem2;
    elem2.e = 0.3;
    elem2.p = 17858.7836;
    elem2.a = elem1.p / (1 - elem1.e * elem1.e);
    elem2.i = 0 * M_PI / 180;
    elem2.flag = 3;
    elem2.w = 30 * M_PI / 180;
    //elem2.W = 90 * M_PI / 180;
    elem2.w_true = 20 * M_PI / 180;
    elem2.mu = 398600.4415;
    auto[delta_v, r, v] = General_plane_change(elem1, elem2);

    COE<double> check = RV2COE(r, v, elem1.mu);
    std::vector<double> h = cross_product(r, v);
    std::vector<double> n = cross_product(std::vector<double>{0, 0, 1}, h);

    std::cout << "norm(n) " << norm(n) << "\n";
    std::cout << h << "\n";
    std::cout << r << "\n" << v << "\n";
    ASSERT_NEAR(check.i, elem2.i, 0.05);
    std::cout << check.i << "\n";
    //ASSERT_NEAR(check.w_true, elem2.w_true, 1e-2);
}



int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
