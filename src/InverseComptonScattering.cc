#include "livpropa/InverseComptonScattering.h"


namespace livpropa {


InverseComptonScattering::InverseComptonScattering(ref_ptr<PhotonField> field, KinematicsMap kinematics, bool havePhotons, double thinning, double limit, unsigned int numberOfSubsteps) {
	setKinematics(kinematics);
	setPhotonField(field);
	setHavePhotons(havePhotons);
	setLimit(limit);
	setThinning(thinning);
	setNumberOfSubsteps(numberOfSubsteps);
	setInteractionTag("EMIC");
}

void InverseComptonScattering::setPhotonField(ref_ptr<PhotonField> field) {
	photonField = field;
	string photonBgName = field->getFieldName();

	string dataPath = "InverseComptonScattering/";
	dataPath += kinematics.getIdentifierForParticle(-11);
	dataPath += "-";
	dataPath += kinematics.getIdentifierForParticle( 11);
	dataPath += "-";
	dataPath += kinematics.getIdentifierForParticle( 22);
	dataPath += "/";
	
	setDescription("InverseComptonScattering: " + photonBgName);
	initRate(getDataPath(dataPath) + "rate_" + photonBgName + ".txt");
	initCumulativeRate(getDataPath(dataPath) + "cdf_" + photonBgName + ".txt");
}

void InverseComptonScattering::setKinematics(KinematicsMap kin) {
	kinematics = kin;
}

void InverseComptonScattering::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

void InverseComptonScattering::setLimit(double limit) {
	this->limit = limit;
}

void InverseComptonScattering::setThinning(double thinning) {
	this->thinning = thinning;
}

void InverseComptonScattering::setInteractionTag(string tag) {
	interactionTag = tag;
}

void InverseComptonScattering::setNumberOfSubsteps(unsigned int n) {
	nSubsteps = n;
}

string InverseComptonScattering::getInteractionTag() const {
	return interactionTag;
}

void InverseComptonScattering::initRate(string filename) {
	std::ifstream infile(filename.c_str());
	if (! infile.good())
		throw runtime_error("InverseComptonScattering: could not open file " + filename);

	// clear previously loaded tables
	tabEnergy.clear();
	tabRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabEnergy.push_back(pow(10, a) * eV);
				tabRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize> ::max(), '\n');
	}
	infile.close();
}

void InverseComptonScattering::initCumulativeRate(string filename) {
	std::ifstream infile(filename.c_str());
	if (! infile.good())
		throw runtime_error("InverseComptonScattering: could not open file " + filename);

	// clear previously loaded tables
	tabE.clear();
	tabs.clear();
	tabCDF.clear();
	
	// skip header
	while (infile.peek() == '#')
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	// read s values in first line
	double a;
	infile >> a; // skip first value
	while (infile.good() and (infile.peek() != '\n')) {
		infile >> a;
		tabs.push_back(pow(10, a) * eV * eV);
	}

	// read all following lines: E, cdf values
	while (infile.good()) {
		infile >> a;
		if (!infile)
			break;  // end of file
		tabE.push_back(pow(10, a) * eV);
		vector<double> cdf;
		for (int i = 0; i < tabs.size(); i++) {
			infile >> a;
			cdf.push_back(a / Mpc);
		}
		tabCDF.push_back(cdf);
	}
	infile.close();
}

// Class to calculate the energy distribution of the ICS photon and to sample from it
class ICSSecondariesEnergyDistribution {
	private:
		vector<vector<double>> data;
		vector<double> s_values;
		size_t Ns;
		size_t Nrer;
		double s_min;
		double s_max;
		double dls;

	public:
		// differential cross-section, see Lee '96 (arXiv:9604098), eq. 23 for x = Ee'/Ee
		double dSigmadE(double x, double beta) {
			double q = ((1 - beta) / beta) * (1 - 1. / x);
			return ((1 + beta) / beta) * (x + 1. / x + 2 * q + q * q);
		}

		// create the cumulative energy distribution of the up-scattered photon
		ICSSecondariesEnergyDistribution() {
			Ns = 1000;
			Nrer = 1000;
			s_min = mec2 * mec2;
			s_max = 1e23 * eV * eV;
			dls = (log(s_max) - log(s_min)) / Ns;
			data = vector<vector<double>>(1000, vector<double>(1000));
			vector<double> data_i(1000);

			// tabulate s bin borders
			s_values = vector<double>(1001);
			for (size_t i = 0; i < Ns + 1; ++i)
				s_values[i] = s_min * exp(i * dls);


			// for each s tabulate cumulative differential cross section
			for (size_t i = 0; i < Ns; i++) {
				double s = s_min * exp((i + 0.5) * dls);
				double beta = (s - s_min) / (s + s_min);
				double x0 = (1 - beta) / (1 + beta);
				double dlx = -log(x0) / Nrer;

				// cumulative midpoint integration
				data_i[0] = dSigmadE(x0, beta) * expm1(dlx);
				for (size_t j = 1; j < Nrer; j++) {
					double x = x0 * exp((j + 0.5) * dlx);
					double dx = exp((j + 1) * dlx) - exp(j * dlx);
					data_i[j] = dSigmadE(x, beta) * dx;
					data_i[j] += data_i[j - 1];
				}
				data[i] = data_i;
			}
		}

		// draw random energy for the up-scattered photon Ep(Ee, s)
		double sample(double Ee, double s) {
			size_t idx = std::lower_bound(s_values.begin(), s_values.end(), s) - s_values.begin();
			vector<double> s0 = data[idx];
			Random &random = Random::instance();
			size_t j = random.randBin(s0) + 1; // draw random bin (upper bin boundary returned)
			double beta = (s - s_min) / (s + s_min);
			double x0 = (1 - beta) / (1 + beta);
			double dlx = -log(x0) / Nrer;
			double binWidth = x0 * (exp(j * dlx) - exp((j - 1) * dlx));
			double Ep = (x0 * exp((j - 1) * dlx) + binWidth) * Ee;
			return std::min(Ee, Ep); // prevent Ep > Ee from numerical inaccuracies
		}
};

void InverseComptonScattering::performInteraction(Candidate* candidate) const {
	int id = candidate->current.getId();

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	if (E < tabE.front() or E > tabE.back())
		return;

	// compute momentum from energy
	double p = kinematics[id]->computeMomentumFromEnergy(E, id);

	// ignore if negative solutions
	if (p < 0)
		return;

	// possible corrections in thresholds
	double sShift = kinematics[id]->getSymmetryBreakingShift(p, id);

	// sample the value of s
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);
	size_t j = random.randBin(tabCDF[i]);
	double s_kin = pow(10, log10(tabs[j]) + (random.rand() - 0.5) * 0.1);
	double s = s_kin + mec2 * mec2 + sShift;

	// sample electron energy after scattering
	static ICSSecondariesEnergyDistribution distribution;
	double Enew = distribution.sample(E, s);
		

	// up-scattered photon
	double Esecondary = E - Enew;

	// avoid possible problems with superluminal case
	if (Esecondary <= 0.)
		return;

	double f = Enew / E;
	if (havePhotons) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
			candidate->addSecondary(22, Esecondary / (1 + z), pos, w, interactionTag);
		}
	}

	// update the primary particle energy; do this after adding the secondary to correctly set the secondary's parent
	candidate->current.setEnergy(Enew / (1 + z));
}

void InverseComptonScattering::process(Candidate* candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	if (E < tabEnergy.front() or (E > tabEnergy.back()))
		return;

	// interaction rate
	double rate = interpolate(E, tabEnergy, tabRate);
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	Random& random = Random::instance();
	int substepCounter = 0;
	do {
		double randDistance = -log(random.rand()) / rate;

		// check for interaction; if it doesn't occur, limit next step
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
			return;
		}
		performInteraction(candidate);

		// performance improvement for interactions near the threshold
		if (substepCounter >= nSubsteps)
			return;

		// repeat with remaining step
		step -= randDistance;
		substepCounter++;

	} while (step > 0);
}

} // namespace livpropa
