%module(directors = "1", threads = "1", allprotected = "1") livpropa


/*************************************************************************************************/
/**                                    SWIG includes                                            **/
/*************************************************************************************************/

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


/*************************************************************************************************/
/**                                   	CRPropa			                                        **/
/*************************************************************************************************/

/* Headers */
%{
	#include "CRPropa.h"
%}

/* Import CRPropa in wrapper */
%import (module = "crpropa") "crpropa.i"


/*************************************************************************************************/
/**                                    LIVpropa preamble                                        **/
/*************************************************************************************************/

%{
	#include "LIVpropa.h"
	using namespace livpropa;
%}


/*************************************************************************************************/
/**                                		Histograms                               				**/
/*************************************************************************************************/

/* To enable access to abstract base class Histogram1D */
%template(Histogram1DRefPtr) crpropa::ref_ptr<livpropa::Histogram1D>;


/*************************************************************************************************/
/**                                		Samplers	                               				**/
/*************************************************************************************************/

%ignore operator livpropa::Sampler*;
%ignore operator livpropa::SamplerList*;

/* To enable access to abstract base class Sampler */

%implicitconv crpropa::ref_ptr<crpropa::SamplerEvents>;
%feature("director") livpropa::SamplerEvents;
%template(SamplerEventsRefPtr) crpropa::ref_ptr<livpropa::SamplerEvents>;

%implicitconv crpropa::ref_ptr<livpropa::SamplerDistribution>;
%template(SamplerDistributionRefPtr) crpropa::ref_ptr<livpropa::SamplerDistribution>;
%feature("director") livpropa::SamplerDistribution;


/*************************************************************************************************/
/**                                		Kinematics                                				**/
/*************************************************************************************************/

/* Ignore list */
%ignore operator livpropa::AbstractKinematics*;
%ignore operator livpropa::LorentzViolating*;

/* To enable access to abstract base class Kinematics */
%implicitconv crpropa::ref_ptr<livpropa::AbstractKinematics>;
%template(AbstractKinematicsRefPtr) crpropa::ref_ptr<livpropa::AbstractKinematics>;
%feature("director") livpropa::AbstractKinematics;

// // Define typemap for crpropa::ref_ptr<livpropa::AbstractKinematics>
// %typemap(out) crpropa::ref_ptr<livpropa::AbstractKinematics> {
// 	$result = SWIG_NewPointerObj(SWIG_as_voidptr($1.get()), SWIGTYPE_p_livpropa__AbstractKinematics, 0 |  0);
// }


/* Print info for SpecialRelativisticKinematics */
%define __STR_SpecialRelativisticKinematics__(SpecialRelativisticKinematics) 
%feature("python:slot", "tp_str", functype = "reprfunc") livpropa::SpecialRelativisticKinematics::_print();
%extend livpropa::SpecialRelativisticKinematics {
	std::string _print() {
		std::ostringstream out;
		out << *$self;
		return out.str();
	}
}
%enddef
__STR_SpecialRelativisticKinematics__(SpecialRelativisticKinematics);

/* Print info for LorentzViolatingKinematicsMonochromatic */
%define __STR_LorentzViolatingKinematicsMonochromatic__(LorentzViolatingKinematicsMonochromatic) 
%feature("python:slot", "tp_str", functype = "reprfunc") livpropa::LorentzViolatingKinematicsMonochromatic::_print();
%extend livpropa::LorentzViolatingKinematicsMonochromatic {
	std::string _print() {
		std::ostringstream out;
		out << *$self;
		return out.str();
	}
}
%enddef
__STR_LorentzViolatingKinematicsMonochromatic__(LorentzViolatingKinematicsMonochromatic);




/* Prevent problems with homonymous functions in different namespace */
%rename(vcMono0_computeMomentumThreshold) livpropa::vc::monoLIV0::computeThresholdMomentum;
%rename(vcMono1_computeMomentumThreshold) livpropa::vc::monoLIV1::computeThresholdMomentum;
%rename(vcMono2_computeMomentumThreshold) livpropa::vc::monoLIV2::computeThresholdMomentum;


/* make Kinematics subscriptable */
// %feature("python:slot", "tp_str", functype = "reprfunc") bsmpropa::Vector4::_print();
// %feature("python:slot", "sq_length", functype = "lenfunc") crpropa::Vector4::__len__;
%feature("python:slot", "mp_subscript", functype = "binaryfunc") livpropa::Kinematics::__getitem__;
%feature("python:slot", "mp_ass_subscript", functype = "objobjargproc") livpropa::Kinematics::__setitem__;
%extend livpropa::Kinematics {
	const crpropa::ref_ptr<livpropa::AbstractKinematics> __getitem__(int i) {
		return (*($self))[i];
	}
}

/* convert unordered_dict to Python dict */
%typemap(out) std::unordered_map<int, crpropa::ref_ptr<livpropa::AbstractKinematics>> (PyObject* obj) %{
	obj = PyDict_New();
	for (const auto& n : $1) {
		// convert int to Python integer
		PyObject* pId = PyLong_FromLong(n.first); // Convert int to Python integer

		// convert crpropa::ref_ptr<livpropa::AbstractKinematics> to Python object
		// PyObject* kin = SWIG_NewPointerObj(new crpropa::ref_ptr<livpropa::AbstractKinematics>(n.second), SWIGTYPE_p_crpropa__ref_ptrT_livpropa__AbstractKinematics_t, SWIG_POINTER_OWN);
		PyObject* kin = SWIG_NewPointerObj(new crpropa::ref_ptr<livpropa::AbstractKinematics>(n.second));

		PyDict_SetItem(obj, pId, kin);
		Py_XDECREF(pId);
		Py_XDECREF(kin);
	}
	$result = SWIG_Python_AppendOutput($result, obj);
%}



/* Instantiate template to ensure that maps of particles-LIV_coefficients are properly handled */
%template(CoefficientsMap) std::unordered_map<int, double>;
%template(ParticleKinematicsMap) std::unordered_map<int, crpropa::ref_ptr<livpropa::AbstractKinematics>>;
// %template(EmissionSpectraTable) std::unordered_map<int, livpropa::EmissionSpectrum>;


/*************************************************************************************************/
/**                                    LIVpropa headers                                        **/
/*************************************************************************************************/


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



/*************************************************************************************************/
/*************************************************************************************************/

%template(VectorInt) std::vector<int>;
%template(VectorFloat) std::vector<float>;
%template(VectorDouble) std::vector<double>;
%template(VectorString) std::vector<std::string>;

/* Ignore list */
%ignore operator<<;
%ignore operator>>;
%ignore *::operator=;


/* Hide warnings */
#pragma SWIG nowarn=312,325,361,389,401,508,509