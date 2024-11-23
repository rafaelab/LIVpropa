#include "livpropa/Kinematics.h"


namespace livpropa {

///////////////////////////////////////////////////////////////////////////////////////////////////

double Kinematics::computeEnergyFromMomentum(const double& p, const int& id) const {
	return sqrt(computeEnergy2FromMomentum(p, id));
}

bool Kinematics::isLorentzInvariant() const {
	if (getNameTag() == "SR") {
		return true;
	} else if (getNameTag() == "LIVmono0") {
		if (static_cast<const MonochromaticLorentzViolatingKinematics<0>*>(this)->getCoefficient() == 0)
			return true;
	} else if (getNameTag() == "LIVmono1") {
		if (static_cast<const MonochromaticLorentzViolatingKinematics<1>*>(this)->getCoefficient() == 0)
			return true;
	} else if (getNameTag() == "LIVmono2") {
		if (static_cast<const MonochromaticLorentzViolatingKinematics<2>*>(this)->getCoefficient() == 0)
			return true;
	}

	return false;
}

bool Kinematics::isLorentzViolating() const {
	return ! isLorentzInvariant();
}

bool Kinematics::isSpecialRelativistic() const {
	return getNameTag() == "SR";
}

bool Kinematics::isLorentzViolatingMonochromatic() const {
	if (getNameTag().find("LIVmono") != std::string::npos)
		return true;

	return false;
}

const SpecialRelativisticKinematics& Kinematics::toSpecialRelativisticKinematics() const {
	return *static_cast<const SpecialRelativisticKinematics*>(this);
}

const MonochromaticLorentzViolatingKinematics<0>& Kinematics::toMonochromaticLorentzViolatingKinematics0() const {
	return *static_cast<const MonochromaticLorentzViolatingKinematics<0>*>(this);
}

const MonochromaticLorentzViolatingKinematics<1>& Kinematics::toMonochromaticLorentzViolatingKinematics1() const {
	return *static_cast<const MonochromaticLorentzViolatingKinematics<1>*>(this);
}

const MonochromaticLorentzViolatingKinematics<2>& Kinematics::toMonochromaticLorentzViolatingKinematics2() const {
	return *static_cast<const MonochromaticLorentzViolatingKinematics<2>*>(this);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

SpecialRelativisticKinematics::SpecialRelativisticKinematics() {
}

SpecialRelativisticKinematics::~SpecialRelativisticKinematics() {
}

string SpecialRelativisticKinematics::getNameTag() const {
	return "SR";
}

double SpecialRelativisticKinematics::getSymmetryBreakingShift(const double& p, const int& id) const {
	return 0;
}

double SpecialRelativisticKinematics::computeEnergy2FromMomentum(const double& p, const int& id) const {
	double m = particleMasses.at(id);
	return pow_integer<2>(p * c_light) + pow_integer<2>(m * c_squared);
}

double SpecialRelativisticKinematics::computeMomentumFromEnergy(const double& E, const int& id) const {
	double m = particleMasses.at(id);
	return sqrt(pow_integer<2>(E) - pow_integer<2>(m * c_squared)) / c_light;
}

double SpecialRelativisticKinematics::getCoefficient() const {
	return 0;
}

string SpecialRelativisticKinematics::getIdentifier() const {
	return "SR";
}

string SpecialRelativisticKinematics::info() const {
	return "SpecialRelativisticKinematics";
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void LorentzViolatingKinematics::setSymmetryBreaking(SymmetryBreaking treatment) {
	symmetryBreaking = treatment;
}

LorentzViolatingKinematics::SymmetryBreaking LorentzViolatingKinematics::getSymmetryBreaking() const {
	return symmetryBreaking;
}

double LorentzViolatingKinematics::selectFinalMomentum(const vector<double>& ps, const double& E, const int& id) const {
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

double LorentzViolatingKinematics::selectFinalMomentumRandom(const vector<double>& ps) {
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

	return -1;
}

double LorentzViolatingKinematics::selectFinalMomentumSmallest(const vector<double>& ps) {
	return *std::min_element(ps.begin(), ps.end());
}

double LorentzViolatingKinematics::selectFinalMomentumLargest(const vector<double>& ps) {
	return *std::max_element(ps.begin(), ps.end());
}

double LorentzViolatingKinematics::selectFinalMomentumAverage(const vector<double>& ps) {
	double pTot = 0;
	for (auto& p : ps) {
		if (p > 0)
			pTot += p;
	}

	return pTot / ps.size();
}

double LorentzViolatingKinematics::selectFinalMomentumClosest(const vector<double>& ps, const double& E, const int& id) {
	SpecialRelativisticKinematics* sr = new SpecialRelativisticKinematics();
	double p0 = sr->computeMomentumFromEnergy(E, id);
	delete sr;

	vector<double> dps;
	for (auto& p : ps) {
		if (p > 0)
			dps.push_back(p - p0);
	}

	return *std::min_element(dps.begin(), dps.end());
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void AbstractMonochromaticLorentzViolatingKinematics::setOrder(int n) {
	order = n;
}

void AbstractMonochromaticLorentzViolatingKinematics::setCoefficient(double coeff) {
	coefficient = coeff;
}

int AbstractMonochromaticLorentzViolatingKinematics::getOrder() const {
	return order;
}

double AbstractMonochromaticLorentzViolatingKinematics::getCoefficient() const {
	return coefficient;
}

string AbstractMonochromaticLorentzViolatingKinematics::getNameTag() const {
	return "LIVmono" + std::to_string(order);
}

string AbstractMonochromaticLorentzViolatingKinematics::getIdentifier() const {
	char s[64];
	sprintf(s, "LIV%i_chi_%+2.1e", order, coefficient);
	return string(s);
}

double AbstractMonochromaticLorentzViolatingKinematics::getSymmetryBreakingShift(const double& p) const {
	return coefficient * pow(p * c_light / energy_planck, order) * pow_integer<2>(p * c_light);
}

string AbstractMonochromaticLorentzViolatingKinematics::info() const {
	string s = "";
	s += "MonochromaticLorentzViolatingKinematics";
	s += "(n = " + std::to_string(order) + "; chi = " + std::to_string(coefficient) + ")"; 
	return s;
}

double AbstractMonochromaticLorentzViolatingKinematics::computeEnergy2FromMomentum(const double& p, const int& id) const {
	auto it = particleMasses.find(id);
	if (it == particleMasses.end()) {
		KISS_LOG_ERROR << "Cannot retrieve particle with id " << id << ". Cannot compute energy from momentum." << endl;
		return 0;
	} 

	double m = (*it).second;
	double ds = getSymmetryBreakingShift(p);

	return pow_integer<2>(p * c_light) + pow_integer<2>(m * c_squared) + ds;
}


///////////////////////////////////////////////////////////////////////////////////////////////////


template<int N>
MonochromaticLorentzViolatingKinematics<N>::MonochromaticLorentzViolatingKinematics(SymmetryBreaking symmetryBreaking) {
	setOrder(N);
	setCoefficient(0.);
	setSymmetryBreaking(symmetryBreaking);
}

template<int N>
MonochromaticLorentzViolatingKinematics<N>::MonochromaticLorentzViolatingKinematics(double chi, SymmetryBreaking symmetryBreaking) {
	setOrder(N);
	setCoefficient(chi);
	setSymmetryBreaking(symmetryBreaking);
}

template<int N>
MonochromaticLorentzViolatingKinematics<N>::~MonochromaticLorentzViolatingKinematics() {
}

template<int N>
double MonochromaticLorentzViolatingKinematics<N>::computeMomentumFromEnergy(const double& E, const int& id) const {
	// LIV for N>2 is not really tested!

	double m = particleMasses.at(id);
	double chi = this->getCoefficient();
	int order = N;

	// vector that will store the solutions for later selection
	vector<double> ps;

	// get vector of coefficients
	unsigned int nCoeffs = order + 3;
	vector<double> coeffsVec;
	coeffsVec.reserve(nCoeffs);
	for (size_t i = 0; i < nCoeffs; i++) {
		coeffsVec.push_back(0);
	}
	coeffsVec[0] = (- pow_integer<2>(E / c_light) + pow_integer<2>(m * c_light)) / c_squared;
	coeffsVec[2] = c_squared;
	coeffsVec[nCoeffs - 1] = chi / pow(c_light / energy_planck, order + 2) * c_squared;

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

	return this->selectFinalMomentum(ps, E, id);
}

template<>
double MonochromaticLorentzViolatingKinematics<0>::computeMomentumFromEnergy(const double& E, const int& id) const {
	double m = particleMasses.at(id);
	double chi = coefficient;
	
	complex<double> sol = sqrt(pow_integer<2>(E) - pow_integer<2>(m * c_squared)) / c_light / sqrt(1 + chi);
	double solRe = sol.real();
	if (solRe < 0)
		solRe = -1;

	vector<double> ps = {solRe};

	return selectFinalMomentum(ps, E, id);
}

template<>
double MonochromaticLorentzViolatingKinematics<1>::computeMomentumFromEnergy(const double& E, const int& id) const {
	double m = particleMasses.at(id);
	double chi = coefficient;

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
		
	return selectFinalMomentum(ps, E, id);
}

template<>
double MonochromaticLorentzViolatingKinematics<2>::computeMomentumFromEnergy(const double& E, const int& id) const {
	double m = particleMasses.at(id);
	double chi = coefficient;
	complex<double> delta = 4 * chi * (E * E - pow_integer<2>(m * c_squared) + pow_integer<2>(energy_planck));
	complex<double> e = energy_planck * sqrt(delta);

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

	return selectFinalMomentum(ps, E, id);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

KinematicsMap::KinematicsMap() {
}

KinematicsMap::KinematicsMap(vector<int> p, vector<ref_ptr<Kinematics>> kin) {
	if (p.size() != kin.size())
		throw std::length_error("Vector of particles and kinematics must have the same length.");

	for (size_t i = 0; i < p.size(); i++) 
		add(p[i], kin[i]);
}

KinematicsMap::KinematicsMap(vector<int> p, ref_ptr<Kinematics> kin) {
	for (size_t i = 0; i < p.size(); i++) 
		add(p[i], kin);
}

KinematicsMap::KinematicsMap(vector<pair<int, ref_ptr<Kinematics>>> kin) {
	for (size_t i = 0; i < kin.size(); i++) 
		add(kin[i].first, kin[i].second);
}

void KinematicsMap::add(const int& particle, ref_ptr<Kinematics> kin) {
	kinematics[particle] = kin;
}

void KinematicsMap::remove(const int& particle) {
	if (! exists(particle))
		throw runtime_error("Cannot retrieve inexistent particle with id " + std::to_string(particle) + ".");

	kinematics.erase(particle);
}

bool KinematicsMap::isLorentzInvariant() const {
	for (auto& kin : kinematics) {
		if (kin.second->getNameTag() != "SR")
			return false;
	}
	return true;
}

bool KinematicsMap::isLorentzViolating() const {
	return ! isLorentzInvariant();
}

bool KinematicsMap::exists(const int& pId) const {
	vector<int> particles = getParticles();
	return std::find(particles.begin(), particles.end(), pId) != particles.end();
}

vector<int> KinematicsMap::getParticles() const {
	vector<int> particles;
	for (auto& kin : kinematics)
		particles.push_back(kin.first);

	return particles;
}

string KinematicsMap::getIdentifierForParticle(const int& pId, bool showParticleId) const {
	char identifier[128] = "";

	const auto& kin = (kinematics.find(pId))->second;
	string kStr = kin->getIdentifier();

	size_t s = 0;
	if (showParticleId)
		s += sprintf(identifier + s, "Id_%+i-%s", pId, kStr.c_str());
	else
		s += sprintf(identifier + s, "%s", kStr.c_str());

	return std::string(identifier);
}

string KinematicsMap::getIdentifier(const std::vector<int>& particles, bool simplify) const {
	string identifier = "";

	// if particles have exactly the same kinematics, return just one identifier.
	if (simplify) {
		int nParticles = particles.size();
		if (nParticles == 1)
			return getIdentifierForParticle(particles[0], false);

		bool simplifiable = true;
		for (size_t i = 1; i < particles.size(); i++) {
			if (getIdentifierForParticle(particles[i - 1], false) != getIdentifierForParticle(particles[i], false)) {
				simplifiable = false;
				break;
			}
		}

		if (simplifiable)
			return getIdentifierForParticle(particles[0], false);
	}

	for (size_t i = 0; i < particles.size(); i++) {
		int pId = particles[i];
		identifier += getIdentifierForParticle(pId);
		if (i < particles.size() - 1)
			identifier += "-";
	}

	return identifier;
}

string KinematicsMap::info() const {
	string s = "KinematicsMap: \n";
	for (auto& kin : kinematics) {
		s += "  . particle " + std::to_string(kin.first) + " ==> " + kin.second->info() + "\n";
	}
	return s;
}

ParticleKinematicsMap KinematicsMap::getParticleKinematicsMap() const {
	return kinematics;
}

const ref_ptr<Kinematics>& KinematicsMap::find(const int& id, bool showWarningInexistent) const {
	bool declared = exists(id);

	if (declared) {
		ParticleKinematicsMapIterator kinIt = kinematics.find(id);
		if (kinIt != kinematics.end()) {
			const ref_ptr<Kinematics>& kin = kinIt->second;
			return kin;
		}
	}

	if (showWarningInexistent) 
		KISS_LOG_WARNING << "Cannot retrieve inexistent particle with id " << id << "." << "Returning special-relativistic kinematics." << endl;

	ref_ptr<Kinematics> specialRelativity = new SpecialRelativisticKinematics();
	return specialRelativity;	
}

const ref_ptr<Kinematics>& KinematicsMap::operator[](const int& pId) {
	return find(pId);
}

ref_ptr<Kinematics> KinematicsMap::operator[](const int& pId) const {
	return find(pId);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

// 	MonochromaticLorentzViolatingKinematics<MonochromaticLIV<N>>* k = new MonochromaticLorentzViolatingKinematics<MonochromaticLIV<N>>(0.);
// 	return *k;
// }

std::ostream& operator<<(std::ostream& os, const Kinematics& kin) {
	os << kin.info();
	return os;
}

std::ostream& operator<<(std::ostream& os, const KinematicsMap& kin) {
	os << "KinematicsMap: " << endl;

	vector<int> particles = kin.getParticles();
	for (auto& p : particles) 
		os << "  . particle " << p << " ==> " << kin[p]->info() << endl;

	return os;
}



// explicit instantiations for required template parameters
template class MonochromaticLorentzViolatingKinematics<0>;
template class MonochromaticLorentzViolatingKinematics<1>;
template class MonochromaticLorentzViolatingKinematics<2>;



} // namespace livpropa