#ifndef LIVPROPA_INVERSECOMPTONSCATTERING_H
#define LIVPROPA_INVERSECOMPTONSCATTERING_H

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include "livpropa/Common.h"
#include "livpropa/Data.h"
#include "livpropa/Kinematics.h"
#include "livpropa/UnitsAndConstants.h"


namespace livpropa {

/**
 @class InverseComptonScattering
 @brief Inverse Compton scattering of electrons with background photons.

 This module simulates inverse Compton scattering of electrons with background photons for several photon fields.
 The upscattered photons are optionally created as secondary particles (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
*/
class InverseComptonScattering: public Module {
	protected:
		ref_ptr<PhotonField> photonField; 	// target photon field
		KinematicsMap kinematics; //!< the type of kinematics (SR, LIV)
		bool havePhotons; // add secondary photons to simulation
		double limit; // limit the step to a fraction of the mean free path
		double thinning; // factor of the thinning (0: no thinning, 1: maximum thinning)
		unsigned int nSubsteps; // do-while loop in process is performed (limiting next step) only nSubsteps times
		string interactionTag;
		vector<double> tabEnergy;  //!< electron energy in [J]
		vector<double> tabRate;  //!< interaction rate in [1/m]
		vector<double> tabE;  //!< electron energy in [J]
		vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
		vector<vector<double>> tabCDF;  //!< cumulative interaction rate

	public:
		InverseComptonScattering(ref_ptr<PhotonField> photonField, KinematicsMap kinematics, bool havePhotons = false, double thinning = 0, double limit = 0.1, unsigned int numberOfSubsteps = 10000);
		void setPhotonField(ref_ptr<PhotonField> photonField);
		void setKinematics(KinematicsMap kin);
		void setHavePhotons(bool havePhotons);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setNumberOfSubsteps(unsigned int n);
		void setInteractionTag(string tag);
		string getInteractionTag() const;
		void initRate(string filename);
		void initCumulativeRate(string filename);
		void process(Candidate* candidate) const;
		void performInteraction(Candidate* candidate) const;
};
/** @}*/

} // namespace livpropa

#endif // LIVPROPA_INVERSECOMPTONSCATTERING_H
