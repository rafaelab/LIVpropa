%module(directors = "1", threads = "1", allprotected = "1") livpropa


%include "attribute.i"
%include "exception.i"
%include "stdint.i"
%include "std_array.i"
%include "std_basic_string.i"
%include "std_container.i"
%include "std_iostream.i"
%include "std_list.i"
%include "std_map.i"
%include "std_set.i"
%include "std_shared_ptr.i"
%include "std_sstream.i"
%include "std_string.i"
%include "std_unordered_map.i"
%include "std_vector.i"
%include "stl.i"
%include "typemaps.i"


/* SWIG exceptions */
%inline %{
	class RangeError {
	};
	class StopIterator {
	};
%}

/* Ignore list */
%ignore operator livpropa::Kinematics*;
%ignore operator livpropa::Sampler*;
%ignore operator livpropa::SamplerList*;


/* Headers */
%{
	#include "CRPropa.h"
	#include "livpropa/UnitsAndConstants.h"
	#include "livpropa/Common.h"
	#include "livpropa/Histogram.h"
	#include "livpropa/Sampler.h"
	#include "livpropa/Data.h"
	#include "livpropa/Kinematics.h"
	#include "livpropa/InverseComptonScattering.h"
	#include "livpropa/PairProduction.h"
	#include "livpropa/PhotonDecay.h"
	#include "livpropa/VacuumCherenkov.h"
	
	using namespace livpropa;
%}


/* Import CRPropa in wrapper */
%import (module = "crpropa") "crpropa.i"


/* To enable access to abstract base class Kinematics */
%implicitconv crpropa::ref_ptr<livpropa::Kinematics>;
%template(KinematicsRefPtr) crpropa::ref_ptr<livpropa::Kinematics>;
%feature("director") livpropa::Kinematics;

%template(Histogram1DRefPtr) crpropa::ref_ptr<livpropa::Histogram1D>;

/* To enable access to abstract base class Sampler */
%implicitconv crpropa::ref_ptr<crpropa::SamplerEvents>;
%feature("director") livpropa::SamplerEvents;
%template(SamplerEventsRefPtr) crpropa::ref_ptr<livpropa::SamplerEvents>;
%implicitconv crpropa::ref_ptr<livpropa::SamplerDistribution>;
%template(SamplerDistributionRefPtr) crpropa::ref_ptr<livpropa::SamplerDistribution>;
%feature("director") livpropa::SamplerDistribution;


/* Fix bug with some SWIG versions */
%feature("notabstract") livpropa::SpecialRelativity;
%feature("notabstract") livpropa::MonochromaticLIV;

/* Prevent problems with homonymous functions in different namespace */
%rename(vcMono0_computeMomentumThreshold) livpropa::vc::monoLIV0::computeThresholdMomentum;
%rename(vcMono1_computeMomentumThreshold) livpropa::vc::monoLIV1::computeThresholdMomentum;
%rename(vcMono2_computeMomentumThreshold) livpropa::vc::monoLIV2::computeThresholdMomentum;


/* Include plugin parts to generate wrappers  */
%include "livpropa/Common.h"
%include "livpropa/Data.h"
%include "livpropa/UnitsAndConstants.h"
%include "livpropa/Histogram.h"
%include "livpropa/Sampler.h"
%include "livpropa/Kinematics.h"
%include "livpropa/InverseComptonScattering.h"
%include "livpropa/PairProduction.h"
%include "livpropa/PhotonDecay.h"
%include "livpropa/VacuumCherenkov.h"


/* Instantiate template to ensure that maps of particles-LIV_coefficients are properly handled */
%template (CoefficientsMap) std::unordered_map<int, double>;




/* Hide warnings */
#pragma SWIG nowarn=312,325,361,389,401,508,509