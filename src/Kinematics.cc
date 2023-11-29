#include "livpropa/Kinematics.h"


namespace livpropa {


double Kinematics::computeEnergyFromMomentum(const double& p, const int& id) const {
	return sqrt(computeEnergy2FromMomentum(p, id));
}



SpecialRelativity::SpecialRelativity() {
}

SpecialRelativity::~SpecialRelativity() {
}

std::string SpecialRelativity::getShortIdentifier() const {
	return "SR";
}

std::string SpecialRelativity::getLocationData(std::vector<int> particles) const {
	return "";
}

double SpecialRelativity::getSymmetryBreakingShift(const double& p, const int& id) const {
	return 0;
}

double SpecialRelativity::computeEnergy2FromMomentum(const double& p, const int& id) const {
	double m = particleMasses.at(id);
	return pow_integer<2>(p * c_light) + pow_integer<2>(m * c_squared);
}

double SpecialRelativity::computeMomentumFromEnergy(const double& E, const int& id) const {
	double m = particleMasses.at(id);
	return sqrt(pow_integer<2>(E) - pow_integer<2>(m * c_squared)) / c_light;
}


MonochromaticLIV::MonochromaticLIV(SymmetryBreakingTreatment symmetryBreakingTreatment) {
}

MonochromaticLIV::MonochromaticLIV(unsigned int n, SymmetryBreakingTreatment symmetryBreakingTreatment) {
	setOrder(n);
}

MonochromaticLIV::MonochromaticLIV(unsigned int order, std::unordered_map<int, double> coeffs, SymmetryBreakingTreatment symmetryBreakingTreatment) {
	setCoefficients(coeffs);
}

MonochromaticLIV::MonochromaticLIV(unsigned int order, std::vector<int> particles, std::vector<double> coeffs, SymmetryBreakingTreatment symmetryBreakingTreatment) {
	if (particles.size() != coeffs.size())
		throw std::length_error("Vectors provided to initialise MonochromaticLIV must be of same size.");

	for (size_t i = 0; i < particles.size(); i++) {
		addCoefficient(particles[i], coeffs[i]);
	}
}

MonochromaticLIV::~MonochromaticLIV() {
}

void MonochromaticLIV::setOrder(unsigned int n) {
	order = n;
}

void MonochromaticLIV::setSymmetryBreakingTreatment(SymmetryBreakingTreatment treatment) {
	symmetryBreakingStrategy = treatment;
}

void MonochromaticLIV::setCoefficients(std::unordered_map<int, double> coeffs) {
	coefficients = coeffs;
}

void MonochromaticLIV::addCoefficient(int particle, double coeff) {
	coefficients.insert({{particle, coeff}});
}

unsigned int MonochromaticLIV::getOrder() const {
	return order;
}

std::unordered_map<int, double> MonochromaticLIV::getCoefficients() const {
	return coefficients;
}

std::string MonochromaticLIV::getShortIdentifier() const {
	return "LIV";
}

std::string MonochromaticLIV::getLocationData(std::vector<int> particles) const {
	char dir[512] = "";

	size_t s = 0;
	if (std::find(particles.begin(), particles.end(), 11) != particles.end()) {
		s += sprintf(dir + s, "chiEl_%+2.1e", getCoefficientForParticle(11));
	} 
	if (std::find(particles.begin(), particles.end(), 22) != particles.end()) {
		s += sprintf(dir + s, "-chiPh_%+2.1e", getCoefficientForParticle(22));
	}
	s += sprintf(dir + s, "-order_%i", order);


	return std::string(dir);
}

std::vector<int> MonochromaticLIV::getParticles() const {
	std::vector<int> particles;
	return particles;
}

double MonochromaticLIV::getCoefficientForParticle(const int& particle) const {
	CoefficientsIterator it = coefficients.find(particle);
	return (*it).second;
}

double MonochromaticLIV::getSymmetryBreakingShift(const double& p, const int& id) const {
	double chi = getCoefficientForParticle(id);
	return chi * pow(p * c_light / energy_planck, getOrder()) * pow_integer<2>(p * c_light);
}

double MonochromaticLIV::computeEnergy2FromMomentum(const double& p, const int& id) const {
	double m = particleMasses.at(id);
	double ds = this->getSymmetryBreakingShift(p, id);
	return pow_integer<2>(p * c_light) + pow_integer<2>(m * c_squared) + ds;
}

double MonochromaticLIV::computeMomentumFromEnergy(const double& E, const int& id) const {
	unsigned int n = getOrder();
	unsigned int nCoeffs = n + 3;
	double m = particleMasses.at(id);
	double chi = getCoefficientForParticle(id);

	std::vector<double> ps;

	if (n == 0) {
		return sqrt(pow_integer<2>(E) - pow_integer<2>(m * c_squared) - chi) / c_light;
	// } 
	// else if (n == 1) {
	// 	return sqrt(pow_integer<2>(E) - pow_integer<2>(m * c_squared)) / c_light;
	} else {
		// get vector of coefficients
		std::vector<double> coeffsVec;
		coeffsVec.reserve(nCoeffs);
		for (size_t i = 0; i < nCoeffs; i++) {
			coeffsVec[i] = 0;
		}
		coeffsVec[0] = - pow_integer<2>(E) + pow_integer<2>(m * c_squared);
		coeffsVec[2] = c_squared;
		coeffsVec[n + 2] = chi / pow(c_light / energy_planck, n + 2) * c_squared;

		// convert std::vector to Eigen vector for root finding
		Eigen::VectorXd coeffs = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(coeffsVec.data(), coeffsVec.size());

		// declare solver and find roots
		Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
		solver.compute(coeffs);
		const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &roots = solver.roots();

		for (auto& root : roots) {
			if (root.real() > 0.)
				ps.push_back(root.real());
		}
	}
	unsigned int nSolutions = ps.size();

	// strategies for treating multiple acceptable solutions
	if (symmetryBreakingStrategy == SymmetryBreakingTreatment::Random) { 
		if (nSolutions == 0) {
			std::cerr << "No solutions found. Momentum cannot be computed from energy in this case." << std::endl;
			return 0;
		} else if (nSolutions == 1) {
			return ps[0];
		} else if (nSolutions == 2) {
			crpropa::Random& random = crpropa::Random::instance();
			if (random.rand() < 0.5)
				return ps[0];
			else
				return ps[1];
		} else {
			crpropa::Random& random = crpropa::Random::instance();

			std::vector<double> edges;
			for (size_t i = 0; i < nSolutions; i++) {
				edges.push_back(1. / nSolutions);
			}
			std::vector<double>::iterator idxMin = std::lower_bound(edges.begin(), edges.end(), random.rand());
			
			return ps[*idxMin];
		}

	} else if (symmetryBreakingStrategy == SymmetryBreakingTreatment::Smallest) {
		std::sort(ps.begin(), ps.end()); 
		return ps[0];
	} else if (symmetryBreakingStrategy == SymmetryBreakingTreatment::Largest) {
		std::sort(ps.begin(), ps.end()); 
		return ps[nSolutions - 1];
	} else if (symmetryBreakingStrategy == SymmetryBreakingTreatment::Average) {
		double pTot = 0;
		for (auto& p : ps)
			pTot += p;
		return pTot / nSolutions;
	} else if (symmetryBreakingStrategy == SymmetryBreakingTreatment::Closest) {
		SpecialRelativity* sr  = new SpecialRelativity();
		double p0 = sr->computeMomentumFromEnergy(E, id);

		std::vector<double> dps;
		for (auto& p : ps) {
			dps.push_back(p - p0);
		}
		std::vector<double>::iterator idx = std::min_element(dps.begin(), dps.end());

		return dps[*idx];
	}


}



} // namespace livpropa