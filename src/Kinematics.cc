#include "livpropa/Kinematics.h"


namespace livpropa {

///////////////////////////////////////////////////////////////////////////////////////////////////

double AbstractKinematics::computeEnergyFromMomentum(const double& p, const int& id) const {
	return sqrt(computeEnergy2FromMomentum(p, id));
}


///////////////////////////////////////////////////////////////////////////////////////////////////

SpecialRelativity::SpecialRelativity() {
}

SpecialRelativity::~SpecialRelativity() {
}

string SpecialRelativity::getShortIdentifier() const {
	return "SR";
}

string SpecialRelativity::getLocationData(const vector<int>& particles) const {
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



///////////////////////////////////////////////////////////////////////////////////////////////////

void LorentzViolating::setSymmetryBreaking(SymmetryBreaking treatment) {
	symmetryBreaking = treatment;
}

LorentzViolating::SymmetryBreaking LorentzViolating::getSymmetryBreaking() const {
	return symmetryBreaking;
}

double LorentzViolating::selectFinalMomentum(const vector<double>& ps, const double& E, const int& id) const {
	unsigned int nSolutions = ps.size();	
	if (nSolutions <= 0) {
		cerr << "No solutions found. Momentum cannot be computed from energy in this case." << endl;
		return 0;
	}

	switch (symmetryBreaking) {
		case SymmetryBreaking::Random:
			return selectFinalMomentumRandom(ps);
		case SymmetryBreaking::Smallest:
			return selectFinalMomentumSmallest(ps);
		case SymmetryBreaking::Largest:
			return selectFinalMomentumLargest(ps);
		case SymmetryBreaking::Average:
			return selectFinalMomentumAverage(ps);
		case SymmetryBreaking::Closest:
			return selectFinalMomentumClosest(ps, E, id);
		default:
			return 0;
	}
}

double LorentzViolating::selectFinalMomentumRandom(const vector<double>& ps) {
	int nSolutions = ps.size();

	if (nSolutions == 1) {
		return ps[0];

	} else if (nSolutions == 2) {
		crpropa::Random& random = crpropa::Random::instance();
		return (random.rand() < 0.5) ? ps[0] : ps[1];

	} else {
		crpropa::Random& random = crpropa::Random::instance();
		vector<double> edges(nSolutions, (double) 1. / nSolutions);
		vector<double>::iterator idxMin = std::lower_bound(edges.begin(), edges.end(), random.rand());
		return ps[*idxMin];

	}

}

double LorentzViolating::selectFinalMomentumSmallest(const vector<double>& ps) {
	return *std::min_element(ps.begin(), ps.end());
}

double LorentzViolating::selectFinalMomentumLargest(const vector<double>& ps) {
	return *std::max_element(ps.begin(), ps.end());
}

double LorentzViolating::selectFinalMomentumAverage(const vector<double>& ps) {
	double pTot = 0;
	for (auto& p : ps)
		pTot += p;

	return pTot / ps.size();
}

double LorentzViolating::selectFinalMomentumClosest(const vector<double>& ps, const double& E, const int& id) {
	SpecialRelativity* sr = new SpecialRelativity();
	double p0 = sr->computeMomentumFromEnergy(E, id);
	
	vector<double> dps;
	for (auto& p : ps)
		dps.push_back(p - p0);

	vector<double>::iterator idx = std::min_element(dps.begin(), dps.end());

	return dps[*idx];
}



///////////////////////////////////////////////////////////////////////////////////////////////////

LorentzViolatingMonochromatic::LorentzViolatingMonochromatic(SymmetryBreaking symmetryBreaking) {
	setSymmetryBreaking(symmetryBreaking);
}

LorentzViolatingMonochromatic::LorentzViolatingMonochromatic(unsigned int n, SymmetryBreaking symmetryBreaking) {
	setOrder(n);
	setSymmetryBreaking(symmetryBreaking);
}

LorentzViolatingMonochromatic::LorentzViolatingMonochromatic(unsigned int n, double chi, SymmetryBreaking symmetryBreaking) {
	setOrder(n);
	setCoefficient(chi);
	setSymmetryBreaking(symmetryBreaking);
}

LorentzViolatingMonochromatic::~LorentzViolatingMonochromatic() {
}

void LorentzViolatingMonochromatic::setOrder(unsigned int n) {
	order = n;
}

void LorentzViolatingMonochromatic::setCoefficient(double coeff) {
	coefficient = coeff;
}

unsigned int LorentzViolatingMonochromatic::getOrder() const {
	return order;
}

double LorentzViolatingMonochromatic::getCoefficient() const {
	return coefficient;
}

string LorentzViolatingMonochromatic::getShortIdentifier() const {
	return "LIV_" + std::to_string(order);
}

double LorentzViolatingMonochromatic::getSymmetryBreakingShift(const double& p) const {
	return coefficient * pow(p * c_light / energy_planck, order) * pow_integer<2>(p * c_light);
}

double LorentzViolatingMonochromatic::computeEnergy2FromMomentum(const double& p, const int& id) const {
	double m = particleMasses.at(id);
	double ds = this->getSymmetryBreakingShift(p);
	return pow_integer<2>(p * c_light) + pow_integer<2>(m * c_squared) + ds;
}

double LorentzViolatingMonochromatic::computeMomentumFromEnergy(const double& E, const int& id) const {
	unsigned int n = order;
	double m = particleMasses.at(id);
	double chi = coefficient;

	// vector that will store the solutions for later selection
	vector<double> ps;
	switch (n) {
		case 0:
			ps = computeMomentumFromEnergy0(E, m, chi);
			break;
		case 1:
			ps = computeMomentumFromEnergy1(E, m, chi);
			break;
		case 2:
			ps = computeMomentumFromEnergy2(E, m, chi);
			break;
		default:
			ps = computeMomentumFromEnergyN(E, m, chi, n);
			break;
	}

	unsigned int nSolutions = ps.size();
	if (nSolutions < 1) {
		cerr << "No solutions could be found. Momentum cannot be computed from energy in this case." << endl;
		return -1;
	}
	// select final momentum according to symmetry breaking treatment
	double p = selectFinalMomentum(ps, E, id);

	return p;
}

vector<double> LorentzViolatingMonochromatic::computeMomentumFromEnergy0(const double& E, const double& m, const double& chi) {
	complex<double> sol = sqrt(pow_integer<2>(E) - pow_integer<2>(m * c_squared)) / c_light / sqrt(1 + chi);
	
	double solRe = sol.real();
	if (solRe < 0)
		solRe = -1;

	vector<double>({solRe});
}

vector<double> LorentzViolatingMonochromatic::computeMomentumFromEnergy1(const double& E, const double& m, const double& chi) {
	complex<double> a = -27 * pow_integer<2>(chi * E) * energy_planck + 2 * pow_integer<3>(energy_planck) + 27 * pow_integer<2>(chi) *  energy_planck * pow_integer<2>(m * c_squared) + sqrt(-4 * pow_integer<6>(energy_planck) + (-27 * pow_integer<2>(chi * E) * energy_planck + 2 * pow_integer<3>(energy_planck) + 27 * pow_integer<2>(chi * m * c_squared) * energy_planck));
	a = pow(a, 1. / 3.);
		
	complex<double> sol1 = - (energy_planck / (3 * chi)) * (1. + pow(2, 1. / 3.) / a) * energy_planck / a - a / 3. / pow(2, 1. / 3.) / chi;
	complex<double> sol2 = - (energy_planck / (3 * chi)) * (1. - energy_planck / pow(2., 2. / 3.) / a) + a / 6. / pow(2., 1. / 3.) / chi;
	
	// third solution is equal to the second (differs in imaginary part), but is considered for averaging algorithm
	complex<double> sol3 = sol2;

	// convert units
	sol1 /= c_light;
	sol2 /= c_light;
	sol3 /= c_light;

	double sol1Re = sol1.real();
	double sol2Re = sol2.real();
	double sol3Re = sol3.real();

	vector<double> ps;
	if (sol1Re > 0)
		ps.push_back(sol1Re);
	if (sol2Re > 0)
		ps.push_back(sol2Re);
	if (sol3Re > 0)
		ps.push_back(sol3Re);
		
	return ps;
}

vector<double> LorentzViolatingMonochromatic::computeMomentumFromEnergy2(const double& E, const double& m, const double& chi) {
	complex<double> e = energy_planck * sqrt(4 * chi * (E * E - pow_integer<2>(m * c_squared) + pow_integer<2>(energy_planck)));
		
	// only two out of the fours solutions lead to positive momenta
	complex<double> sol1 = sqrt((- pow_integer<2>(energy_planck) - energy_planck * e) / chi / 2.);
	complex<double> sol2 = sqrt((- pow_integer<2>(energy_planck) + energy_planck * e) / chi / 2.);

	double sol1Re = sol1.real();
	double sol2Re = sol2.real();
	
	vector<double> ps;
	if (sol1Re > 0)
		ps.push_back(sol1Re);
	if (sol2Re > 0)
		ps.push_back(sol2Re);

	return ps;
}

vector<double> LorentzViolatingMonochromatic::computeMomentumFromEnergyN(const double& E, const double& m, const double& chi, const unsigned int& n) {
	// UNTESTED!!!!

	// get vector of coefficients
	unsigned int nCoeffs = n + 3;
	vector<double> coeffsVec;
	coeffsVec.reserve(nCoeffs);
	for (size_t i = 0; i < nCoeffs; i++) {
		coeffsVec.push_back(0);
	}
	coeffsVec[0] = (- pow_integer<2>(E / c_light) + pow_integer<2>(m * c_light)) / c_squared;
	coeffsVec[2] = c_squared;
	coeffsVec[nCoeffs - 1] = chi / pow(c_light / energy_planck, n + 2) * c_squared;

	// convert vector to Eigen vector for root finding
	Eigen::VectorXd coeffs = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(coeffsVec.data(), coeffsVec.size());

	// declare solver and find roots
	Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
	solver.compute(coeffs);
	const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &roots = solver.roots();

	vector<double> ps;
	for (auto& root : roots) {
		if (root.real() > 0.)
			ps.push_back(root.real());
	}

	return ps;
}



///////////////////////////////////////////////////////////////////////////////////////////////////

MonochromaticLIV::MonochromaticLIV(SymmetryBreaking symmetryBreaking) {
}

MonochromaticLIV::MonochromaticLIV(unsigned int n, SymmetryBreaking symmetryBreaking) {
	setOrder(n);
}

MonochromaticLIV::MonochromaticLIV(unsigned int n, double chi, SymmetryBreaking symmetryBreaking) {
	setOrder(n);
	addCoefficient(0, chi);
}

MonochromaticLIV::MonochromaticLIV(unsigned int order, unordered_map<int, double> coeffs, SymmetryBreaking symmetryBreaking) {
	setCoefficients(coeffs);
}

MonochromaticLIV::MonochromaticLIV(unsigned int order, vector<int> particles, vector<double> coeffs, SymmetryBreaking symmetryBreaking) {
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

void MonochromaticLIV::setCoefficients(unordered_map<int, double> coeffs) {
	coefficients = coeffs;
}

void MonochromaticLIV::addCoefficient(int particle, double coeff) {
	coefficients.insert({{particle, coeff}});
}

unsigned int MonochromaticLIV::getOrder() const {
	return order;
}

unordered_map<int, double> MonochromaticLIV::getCoefficients() const {
	return coefficients;
}

string MonochromaticLIV::getShortIdentifier() const {
	return "LIV";
}

string MonochromaticLIV::getLocationData(const vector<int>& particles) const {
	char dir[512] = "";

	size_t s = 0;
	if (std::find(particles.begin(), particles.end(), 11) != particles.end()) {
		s += sprintf(dir + s, "chiEl_%+2.1e", getCoefficientForParticle(11));
	} 
	if (std::find(particles.begin(), particles.end(), 22) != particles.end()) {
		s += sprintf(dir + s, "-chiPh_%+2.1e", getCoefficientForParticle(22));
	}
	s += sprintf(dir + s, "-order_%i", order);


	return string(dir);
}

vector<int> MonochromaticLIV::getParticles() const {
	vector<int> particles;
	return particles;
}

double MonochromaticLIV::getCoefficientForParticle(const int& particle) const {
	CoefficientsIterator it = coefficients.find(particle);

	if (it == coefficients.end()) {
		CoefficientsIterator it0 = coefficients.find(0);
		if (it0 == coefficients.end()) {
			cout << "Particle not found in list. Returning chi=0." << endl;
			return 0;
		} else {
			return (*it0).second;
		}
	}
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

	vector<double> ps;

	if (n == 0) {
		complex<double> sol = sqrt(pow_integer<2>(E) - pow_integer<2>(m * c_squared)) / c_light / sqrt(1 + chi);
		double solRe = sol.real();
		if (solRe < 0)
			return -1;
		return solRe;
		
	} else if (n == 1) {

		complex<double> a = -27 * pow_integer<2>(chi * E) * energy_planck + 2 * pow_integer<3>(energy_planck) + 27 * pow_integer<2>(chi) *  energy_planck * pow_integer<2>(m * c_squared) + sqrt(-4 * pow_integer<6>(energy_planck) + (-27 * pow_integer<2>(chi * E) * energy_planck + 2 * pow_integer<3>(energy_planck) + 27 * pow_integer<2>(chi * m * c_squared) * energy_planck));
		a = pow(a, 1. / 3.);
		
		complex<double> sol1 = - (energy_planck / (3 * chi)) * (1. + pow(2, 1. / 3.) / a) * energy_planck / a - a / 3. / pow(2, 1. / 3.) / chi;
		complex<double> sol2 = - (energy_planck / (3 * chi)) * (1. - energy_planck / pow(2., 2. / 3.) / a) + a / 6. / pow(2., 1. / 3.) / chi;
		
		// third solution is equal to the second (differs in imaginary part), but is considered for averaging algorithm
		complex<double> sol3 = sol2;

		// convert units
		sol1 /= c_light;
		sol2 /= c_light;
		sol3 /= c_light;

		double sol1Re = sol1.real();
		double sol2Re = sol2.real();
		double sol3Re = sol3.real();

		if (sol1Re > 0)
			ps.push_back(sol1Re);
		if (sol2Re > 0)
			ps.push_back(sol2Re);
		if (sol3Re > 0)
			ps.push_back(sol3Re);

	} else if (n == 2) {
		complex<double> e = energy_planck * sqrt(4 * chi * (E * E - pow_integer<2>(m * c_squared) + pow_integer<2>(energy_planck)));
		
		// only two out of the fours solutions lead to positive momenta
		complex<double> sol1 = sqrt((- pow_integer<2>(energy_planck) - energy_planck * e) / chi / 2.);
		complex<double> sol2 = sqrt((- pow_integer<2>(energy_planck) + energy_planck * e) / chi / 2.);

		double sol1Re = sol1.real();
		double sol2Re = sol2.real();
		
		if (sol1Re > 0)
			ps.push_back(sol1Re);
		if (sol2Re > 0)
			ps.push_back(sol2Re);

	} else {
		// get vector of coefficients
		vector<double> coeffsVec;
		coeffsVec.reserve(nCoeffs);
		for (size_t i = 0; i < nCoeffs; i++) {
			coeffsVec.push_back(0);
		}
		coeffsVec[0] = (- pow_integer<2>(E / c_light) + pow_integer<2>(m * c_light)) / c_squared;
		coeffsVec[2] = c_squared;
		coeffsVec[nCoeffs - 1] = chi / pow(c_light / energy_planck, n + 2) * c_squared;

		// convert vector to Eigen vector for root finding
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
	if (nSolutions < 1)
		return -1;

	// strategies for treating multiple acceptable solutions
	if (symmetryBreaking == SymmetryBreaking::Random) { 
		if (nSolutions == 0) {
			cerr << "No solutions found. Momentum cannot be computed from energy in this case." << endl;
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

			vector<double> edges;
			for (size_t i = 0; i < nSolutions; i++) {
				edges.push_back(1. / nSolutions);
			}
			vector<double>::iterator idxMin = std::lower_bound(edges.begin(), edges.end(), random.rand());
			
			return ps[*idxMin];
		}

	} else if (symmetryBreaking == SymmetryBreaking::Smallest) {
		std::sort(ps.begin(), ps.end()); 
		return ps[0];
	} else if (symmetryBreaking == SymmetryBreaking::Largest) {
		std::sort(ps.begin(), ps.end()); 
		return ps[nSolutions - 1];
	} else if (symmetryBreaking == SymmetryBreaking::Average) {
		double pTot = 0;
		for (auto& p : ps)
			pTot += p;
		return pTot / nSolutions;
	} else if (symmetryBreaking == SymmetryBreaking::Closest) {
		SpecialRelativity* sr  = new SpecialRelativity();
		double p0 = sr->computeMomentumFromEnergy(E, id);

		vector<double> dps;
		for (auto& p : ps) {
			dps.push_back(p - p0);
		}
		vector<double>::iterator idx = std::min_element(dps.begin(), dps.end());

		return dps[*idx];
	}


}



} // namespace livpropa