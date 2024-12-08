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
%include "std_pair.i"
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


/* C++ containers */
%template() std::pair<double, double>;
%template(PairVector) std::vector<std::pair<double, double>>;


/*************************************************************************************************/
/**	                        			 NumPy interface  		                                **/
/*************************************************************************************************/

%template(VectorInt) std::vector<int>;
%template(VectorFloat) std::vector<float>;
%template(VectorDouble) std::vector<double>;
%template(VectorString) std::vector<std::string>;

%feature("director:except") {
	if ($error != NULL) {
		PyObject* ptype;
		PyObject* pvalue; 
		PyObject* ptraceback;
		PyErr_Fetch(&ptype, &pvalue, &ptraceback);
		PyErr_Restore(ptype, pvalue, ptraceback);
		PyErr_Print();
		Py_Exit(1);
	}
}

/* Exceptions for Python lists and iterators */
%exception __next__ {
	try {
		$action
	} catch (StopIterator) {
		PyErr_SetString(PyExc_StopIteration, "End of iterator");
		return NULL;
	}
}

%exception __getitem__ {
	try {
		$action
	} catch (RangeError) {
		SWIG_exception(SWIG_IndexError, "Index out of bounds");
		return NULL;
	}
};

%{
	#define SWIG_FILE_WITH_INIT
	#include <vector>
	#include <numpy/arrayobject.h>
	#include "numpy/ufuncobject.h"
%}
%include "numpy.i"
%init %{
	import_array();
	import_ufunc();
%}

// typemap for converting std::vector<double> to numpy array
%typemap(out) const std::vector<double>& {
	npy_intp size = $1.size();

	PyObject* obj = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
	if (! obj) {
		SWIG_exception_fail(SWIG_RuntimeError, "Unable to create numpy array");
	}

	double* data = static_cast<double*>(PyArray_DATA((PyArrayObject*) obj));
	for (npy_intp i = 0; i < size; ++i) {
		data[i] = $1[i];
	}

	$result = obj;
}

// typemap for converting numpy array to std::vector<double>
%typemap(in) (std::vector<double>& vec) {
	PyArrayObject* array = (PyArrayObject*) PyArray_FROM_OTF($input, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	if (! array) {
		SWIG_exception_fail(SWIG_TypeError, "Expected a numpy array of type float64");
	}

	npy_intp size = PyArray_SIZE(array);
	double* data = static_cast<double*>(PyArray_DATA(array));
	vec.assign(data, data + size)
	Py_DECREF(array);

	$1 = vec;
}

%apply(double* INPLACE_ARRAY1, int DIM1) { 
	(double *c, int len_c) 
};

%apply(double* ARGOUT_ARRAY1, int DIM1) {
	(double* rangevec, int n)
};

%apply unsigned int &OUTPUT {
	unsigned int&
}; 

/* Python slots */
%feature("python:slot", "sq_length", functype = "lenfunc") __len__;
%feature("python:slot", "mp_subscript", functype = "binaryfunc") __getitem__;
%feature("python:slot", "tp_iter", functype = "unaryfunc") __iter__;
%feature("python:slot", "tp_iternext", functype = "iternextfunc") __next__;

%typemap(directorin, numinputs = 1) (const double* v) {
	npy_intp dim = v.size();
	$input = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, (void*) $1);
}

%typemap(directorin, numinputs = 1) (const std::vector<double>& v) {
	npy_intp dim = v.size();
	$input = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, (void*) $1);
}

%fragment("NumPy_Fragments");
%fragment("NumPy_Macros");
%numpy_typemaps(double, NPY_DOUBLE, size_t)
%numpy_typemaps(int, NPY_INT, size_t)

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

/* Include plugin parts to generate wrappers  */
%include "livpropa/Common.h"
%include "livpropa/Data.h"
%include "livpropa/UnitsAndConstants.h"
%include "livpropa/Histogram.h"
%include "livpropa/Sampler.h"
%include "livpropa/Weighter.h"
%include "livpropa/Kinematics.h"
%include "livpropa/InverseComptonScattering.h"
%include "livpropa/PairProduction.h"
%include "livpropa/PhotonDecay.h"
%include "livpropa/VacuumCherenkov.h"


/*************************************************************************************************/
/**                                		Histograms                               				**/
/*************************************************************************************************/

/* Ignore list */
%ignore operator livpropa::Bin1D*;
%ignore operator livpropa::Histogram1D*;

/* Rename automatically generate enum class */
%rename(LogBase_e) livpropa::LogBase_e;
%rename(LogBase_2) livpropa::LogBase_two;
%rename(LogBase_10) livpropa::LogBase_ten;

/* To enable access to abstract base class Bin1D */
%implicitconv crpropa::ref_ptr<livpropa::Bin1D>;
%template(Bin1DRefPtr) crpropa::ref_ptr<livpropa::Bin1D>;
%feature("director") livpropa::Bin1D;

/* To enable access to abstract base class Histogram1D */
%implicitconv crpropa::ref_ptr<livpropa::Histogram1D>;
%template(Histogram1DRefPtr) crpropa::ref_ptr<livpropa::Histogram1D>;
%feature("director") livpropa::Histogram1D;


/*************************************************************************************************/
/**	                          Samplers & Weighters  			                                    **/
/*************************************************************************************************/

/* Ignore list */
%ignore operator livpropa::Sampler*;
%ignore operator livpropa::Weighter*;
%ignore operator livpropa::WeighterList*;

/* To enable access to abstract base class  */
%implicitconv crpropa::ref_ptr<livpropa::Sampler>;
%template(SamplerRefPtr) crpropa::ref_ptr<livpropa::Sampler>;
%feature("director") livpropa::Sampler;

/* allow mutables */
%feature("immutable", "1") livpropa::Nested;
%feature("immutable", "0") livpropa::NestedSampler::livePoints;
%feature("immutable", "0") livpropa::NestedSampler::liveLikelihoods;
%feature("immutable", "0") livpropa::NestedSampler::logEvidence;
%feature("immutable", "0") livpropa::NestedSampler::logWeight;

%feature("immutable", "1") livpropa::MCMCSampler;
%feature("immutable", "0") livpropa::MCMCSampler::currentSample;
%feature("immutable", "0") livpropa::MCMCSampler::currentWeight;

%feature("immutable", "1") livpropa::AdaptiveMCMCSampler;
%feature("immutable", "0") livpropa::AdaptiveMCMCSampler::stepSize;
%feature("immutable", "0") livpropa::AdaptiveMCMCSampler::currentSample;
%feature("immutable", "0") livpropa::AdaptiveMCMCSampler::currentWeight;
%feature("immutable", "0") livpropa::AdaptiveMCMCSampler::acceptedSamples;
%feature("immutable", "0") livpropa::AdaptiveMCMCSampler::acceptanceRate;


/* To enable access to abstract base class  */
%implicitconv crpropa::ref_ptr<livpropa::Weighter>;
%template(WeighterRefPtr) crpropa::ref_ptr<livpropa::Weighter>;
%feature("director") crpropa::Weighter;


/*************************************************************************************************/
/**                                		Kinematics                                				**/
/*************************************************************************************************/

/* Ignore list */
%ignore operator livpropa::Kinematics*;

/* To enable access to abstract base class Kinematics */
%implicitconv crpropa::ref_ptr<livpropa::Kinematics>;
%template(KinematicsRefPtr) crpropa::ref_ptr<livpropa::Kinematics>;
%feature("director") livpropa::Kinematics;

%template(MonochromaticLorentzViolatingKinematics0) livpropa::MonochromaticLorentzViolatingKinematics<0>;
%template(MonochromaticLorentzViolatingKinematics1) livpropa::MonochromaticLorentzViolatingKinematics<1>;
%template(MonochromaticLorentzViolatingKinematics2) livpropa::MonochromaticLorentzViolatingKinematics<2>;

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
%define __STR_AbstractMonochromaticLorentzViolatingKinematics__(AbstractMonochromaticLorentzViolatingKinematics) 
%feature("python:slot", "tp_str", functype = "reprfunc") livpropa::AbstractMonochromaticLorentzViolatingKinematics::_print();
%extend livpropa::AbstractMonochromaticLorentzViolatingKinematics {
	std::string _print() {
		std::ostringstream out;
		out << *$self;
		return out.str();
	}
}
%enddef
__STR_AbstractMonochromaticLorentzViolatingKinematics__(AbstractMonochromaticLorentzViolatingKinematics);


/* make KinematicsMap subscriptable */
%feature("python:slot", "mp_subscript", functype = "binaryfunc") livpropa::KinematicsMap::__getitem__;
%feature("python:slot", "mp_ass_subscript", functype = "objobjargproc") livpropa::KinematicsMap::__setitem__;
%extend livpropa::KinematicsMap {
	const crpropa::ref_ptr<livpropa::Kinematics> __getitem__(int i) {
		return (*($self))[i];
	}
}

/* print info Kinematics */
%define __STR_Kinematics__(AbstractMonochromaticLorentzViolatingKinematics) 
%feature("python:slot", "tp_str", functype = "reprfunc") livpropa::KinematicsMap::_print();
%extend livpropa::Kinematics {
	std::string _print() {
		std::ostringstream out;
		out << (*$self).info();
		return out.str();
	}
}
%enddef
__STR_Kinematics__(KinematicsMap);


/* convert unordered_dict to Python dict */
%typemap(out) std::unordered_map<int, crpropa::ref_ptr<livpropa::Kinematics>>(PyObject* obj) %{
	obj = PyDict_New();
	for (const auto& n : $1) {
		// convert int to Python integer
		PyObject* pId = PyLong_FromLong(n.first); // Convert int to Python integer

		// convert crpropa::ref_ptr<livpropa::Kinematics> to Python object
		PyObject* kin = SWIG_NewPointerObj(new crpropa::ref_ptr<livpropa::Kinematics>(n.second));

		PyDict_SetItem(obj, pId, kin);
		Py_XDECREF(pId);
		Py_XDECREF(kin);
	}
	$result = SWIG_Python_AppendOutput($result, obj);
%}


/*************************************************************************************************/
/**	                                	Interactions  		                                    **/
/*************************************************************************************************/

/* Rename automatically generate enum class */
%rename(VacuumCherenkovSpectrumDefault) VacuumCherenkovSpectrum_Default;
%rename(VacuumCherenkovSpectrumFull) VacuumCherenkovSpectrum_Full;
%rename(VacuumCherenkovSpectrumStep) VacuumCherenkovSpectrum_Step;
%rename(VacuumCherenkovSpectrumAbsent) VacuumCherenkovSpectrum_Absent;


/*************************************************************************************************/
/*************************************************************************************************/




%clear(double* vector, int length);


/* Ignore list */
%ignore operator<<;
%ignore operator>>;
%ignore *::operator=;


/* Hide warnings */
#pragma SWIG nowarn=302,312,315,325,361,389,401,508,509