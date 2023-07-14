#ifndef LIVPROPA_INVERSECOMPTONSCATTERINGLIV_H
#define LIVPROPA_INVERSECOMPTONSCATTERINGLIV_H

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

#include <crpropa/Candidate.h>
#include <crpropa/Common.h>
#include <crpropa/Module.h>
#include <crpropa/PhotonBackground.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>

#include "livpropa/Data.h"
#include "livpropa/LorentzSymmetry.h"
#include "livpropa/UnitsAndConstants.h"

using crpropa::Candidate;
using crpropa::Module;
using crpropa::PhotonField;
using crpropa::Random;
using crpropa::Vector3d;
using crpropa::ref_ptr;
using crpropa::pow_integer;
using crpropa::interpolate;
using crpropa::closestIndex;

namespace livpropa {

/**
 @class InverseComptonScatteringLIV
 @brief Inverse Compton scattering of electrons with background photons.

 This module simulates inverse Compton scattering of electrons with background photons for several photon fields.
 The upscattered photons are optionally created as secondary particles (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
*/
class InverseComptonScatteringLIV: public Module {
private:
	ref_ptr<PhotonField> photonField;
	bool havePhotons;
	double limit;
	double thinning;
	std::string interactionTag = "EMIC";

	// tabulated interaction rate 1/lambda(E)
	std::vector<double> tabEnergy;  //!< electron energy in [J]
	std::vector<double> tabRate;  //!< interaction rate in [1/m]
	
	// tabulated CDF(s_kin, E) = cumulative differential interaction rate
	std::vector<double> tabE;  //!< electron energy in [J]
	std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
	std::vector< std::vector<double> > tabCDF;  //!< cumulative interaction rate

public:
	/** Constructor
	 @param photonField		target photon field
	 @param havePhotons		if true, add secondary photons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
	 */
	InverseComptonScatteringLIV(ref_ptr<PhotonField> photonField, bool havePhotons = false, double thinning = 0, double limit = 0.1);

	// set the target photon field 
	void setPhotonField(ref_ptr<PhotonField> photonField);

	// decide if secondary photons are added to the simulation
	void setHavePhotons(bool havePhotons);

	/** limit the step to a fraction of the mean free path
	 @param limit	fraction of the mean free path, should be between 0 and 1
	*/
	void setLimit(double limit);

	/** Apply thinning with a given thinning factor
	 * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
	 */
	void setThinning(double thinning);

	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */
	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;

	void initRate(std::string filename);
	void initCumulativeRate(std::string filename);

	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;
};
/** @}*/

} // namespace livpropa

#endif // LIVPROPA_INVERSECOMPTONSCATTERINGLIV_H
