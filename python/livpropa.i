%module(directors = "1", threads = "1", allprotected = "1") livpropa


%include "attribute.i"
%include "exception.i"
%include "stdint.i"
%include "std_array.i"
%include "std_container.i"
%include "std_iostream.i"
%include "std_list.i"
%include "std_map.i"
%include "std_set.i"
%include "std_shared_ptr.i"
%include "std_string.i"
%include "std_vector.i"
%include "stl.i"
%include "typemaps.i"


// /* Ignore list */
// %ignore operator alpinist::PlasmaDensity*;


/* Headers */
%{
	#include "CRPropa.h"
	#include "livpropa/Data.h"
	#include "livpropa/UnitsAndConstants.h"
	#include "livpropa/LorentzSymmetry.h"
	#include "livpropa/PairProductionLIV.h"
	
	using namespace livpropa;
%}


/* Import CRPropa in wrapper */
%import (module = "crpropa") "crpropa.i"


// /* Templates for ref_ptr */
// %implicitconv crpropa::ref_ptr<alpinist::PlasmaDensity>;
// %template(PlasmaDensityRefPtr) crpropa::ref_ptr<alpinist::PlasmaDensity>;
// %feature("director") alpinist::PlasmaDensity;


/* Include plugin parts to generate wrappers for */
%include "livpropa/Data.h"
%include "livpropa/UnitsAndConstants.h"
%include "livpropa/LorentzSymmetry.h"
%include "livpropa/PairProductionLIV.h"



/* Hide warnings */
#pragma SWIG nowarn=312,325,361,389,401,508,509