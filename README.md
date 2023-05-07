Library, that estimates delta-v

Input:
Keplerian elements of initial and final orbits:
* p - semilatus rectum
* a - semimajor axis
* e - eccentricity
* i - inclination
* W - right ascension (for inclined orbits)
* w - argument of perigee (for elliptic inclined orbits)
* nu - true anomaly (for now optional)
* u - argument of latitude (for circular inclined orbits)
* lam_true - true longitude (for circular equatorial orbits)
* w_true - true longitude of periapsis (for elliptic equatorial orbits)
* mu - gravitational parameter
* flag - determines type of orbit:1 - circular equatorial, 2 - circular inclined, 3 - elliptic equatorial, 4 - elliptic inclined

Output: delta-v

Orbital_elements_convertion.h
1) RV2COE:
Function converts RV vectors to Keplerian elements

2) COE2RV:
Function converts Keplerian elemnts to RV vectors

Orbital_maneuvers.h
1) Hohmann_transfer:
Hohmann algorithm to estimate delta-v

2) Bi_elliptic_transfer_circular_orbits:
Bi-elliptic algorithm to estimate delta-v for circular orbits

3) Bi_elliptic_transfer_elliptic_orbits:
Bi-elliptic transfer algorithm to estimate delta-v for elliptic orbits

4) Two_impulse_transfer_elliptic_orbits:
Two impulse transfer algorithm to estimate delta-v for elliptic orbits

5) Inclination_only_transfer:
Inclination only plane change transfer algorithm to estimate delta-v

6) General_plane_change:
General plane change algorithm to estimate delta-v

* Finding planes intersection vectors
* Finding 2 points, that belongs to initial orbit and final plane
* Considering, that V vector only changes its direction(without changing its magnitude), estimating delta-v by calculating angle between planes(difference between inclination angles)
* Returning delta-v and RV vectors in order to find transfer orbit, which plane coincide with final orbit's plane

7) General_transfer:
Combining general plane change transfer with coplanar transfer
