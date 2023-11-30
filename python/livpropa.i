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


%begin %{
	#define SWIG_PYTHON_CAST_MODE
%}

/* Headers */
%{
	#include "CRPropa.h"
	#include "livpropa/UnitsAndConstants.h"
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


/* Fix bug with some SWIG versions */
%feature("notabstract") livpropa::SpecialRelativity;
%feature("notabstract") livpropa::MonochromaticLIV;


/* Include plugin parts to generate wrappers  */
%include "livpropa/Data.h"
%include "livpropa/UnitsAndConstants.h"
%include "livpropa/Kinematics.h"
%include "livpropa/InverseComptonScattering.h"
%include "livpropa/PairProduction.h"
%include "livpropa/PhotonDecay.h"
%include "livpropa/VacuumCherenkov.h"

/* Instantiate template to ensure that maps of particles-LIV_coefficients are propertly handled */
%template (CoefficientsMap) std::unordered_map<int, double>;

/* Hide warnings */
#pragma SWIG nowarn=312,325,361,389,401,508,509