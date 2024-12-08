#ifndef LIVPROPA_COMMON_H
#define LIVPROPA_COMMON_H

#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <crpropa/Candidate.h>
#include <crpropa/Common.h>
#include <crpropa/Module.h>
#include <crpropa/PhotonBackground.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Vector3.h>


namespace livpropa {


// import from std namespace
using std::cerr;
using std::cout;
using std::endl;
using std::complex;
using std::string;
using std::vector;
using std::pair;
using std::unordered_map;
using std::runtime_error;


// import from CRPropa namespace
using crpropa::Candidate;
using crpropa::Module;
using crpropa::PhotonField;
using crpropa::Random;
using crpropa::Referenced;
using crpropa::Vector3d;
using crpropa::ref_ptr;
using crpropa::pow_integer;
using crpropa::interpolate;
using crpropa::interpolateEquidistant;
using crpropa::closestIndex;


/**
 @brief Calculate the logarithm of a number with a given base.
 */
inline double logBase(double x, double b) {
	return log(x) / log(b);
}

/**
  @brief Linear interpolation between two points.
 */
inline double twoPointExtrapolation(const double& xi, const double& x1, const double& y1, const double& x2, const double& y2) {
	if (x1 == x2) {
		cout << "Error: x1 and x2 cannot be the same value. Returning NaN." << endl;
		return std::nan("");
	}

	double dydx = (y2 - y1) / (x2 - x1);
	double yi = y1 + dydx * (xi - x1);

	return yi;
}


} // livpropa


#endif // LIVPROPA_COMMON_H
