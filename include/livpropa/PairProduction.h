#ifndef LIVPROPA_PAIRPRODUCTION_H
#define LIVPROPA_PAIRPRODUCTION_H

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
 @class PairProduction
 @brief Breit-Wheeler pair production of photons with background photons considering Lorentz invariance violation (LIV).

 This module simulates electron-pair production of photons with background photons:
 	gamma + gamma_b --> e+ + e- 
 The resulting electron-positron pair is optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 */
class PairProduction: public Module {
	protected:
		ref_ptr<PhotonField> photonField; 	// target photon field
		Kinematics kinematics; //!< the type of kinematics (SR, LIV)
		bool haveElectrons; // add secondary electrons to simulation
		double limit; // limit the step to a fraction of the mean free path
		double thinning; // factor of the thinning (0: no thinning, 1: maximum thinning)
		string interactionTag;  //!< tag corresponding to this interaction
		vector<double> tabEnergy;  //!< electron energy in [J]
		vector<double> tabRate;  //!< interaction rate in [1/m]
		vector<double> tabE;  //!< electron energy in [J]
		vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
		vector<vector<double>> tabCDF;  //!< cumulative interaction rate

	public:
		PairProduction();
		PairProduction(ref_ptr<PhotonField> photonField, Kinematics kinematics, bool haveElectrons = false, double thinning = 0, double limit = 0.1);
		void setPhotonField(ref_ptr<PhotonField> photonField);
		void setHaveElectrons(bool haveElectrons);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setKinematics(Kinematics kin);
		void setInteractionTag(string tag);
		string getInteractionTag() const;
		void initRate(string filename);
		void initCumulativeRate(string filename);
		void performInteraction(Candidate* candidate) const;
		void process(Candidate* candidate) const;
};


} // namespace livpropa

#endif // LIVPROPA_PAIRPRODUCTION_H
