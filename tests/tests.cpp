#include "gtest/gtest.h"
#include "../src/Orbital_elements_convertion.h"
#include "../src/Orbital_maneuvers.h"

TEST(ORBITAL_MANEUVERS, RV2COE_ELLIPTIC_INCLINED) {
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
    auto[r, v] = COE2RV(elem);
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


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
