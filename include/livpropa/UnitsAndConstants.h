#ifndef LIVPROPA_UNITSANDCONSTANTS_H
#define LIVPROPA_UNITSANDCONSTANTS_H

#include <unordered_map>

#include <crpropa/Units.h>

#include "livpropa/Common.h"


namespace livpropa {

// constants
using crpropa::mass_electron;
using crpropa::c_light;
using crpropa::c_squared;
using crpropa::eplus;
using crpropa::h_planck;
using crpropa::alpha_finestructure;

// units
using crpropa::eV;
using crpropa::kpc;
using crpropa::Mpc;


// LIVpropa-specific definitions
static const double h_dirac = h_planck / 2. / M_PI;
static const double mass_planck = 2.176434e-8; // in kg
static const double energy_planck = mass_planck * c_squared;

static const double mec2 = mass_electron * c_squared;

// particle masses
static const unordered_map<int, double> particleMasses = {
	{ 11, mass_electron}, 
	{-11, mass_electron},
	{ 22, 0.}
};

} // namespace livpropa


#endif // LIVPROPA_UNITSANDCONSTANTS_H