#ifndef LIVPROPA_COMMON_H
#define LIVPROPA_COMMON_H

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
using crpropa::closestIndex;



inline double logBase(double x, double b) {
    return log(x) / log(b);
}



} // livpropa


#endif // LIVPROPA_COMMON_H
