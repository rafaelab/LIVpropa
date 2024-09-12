#include "livpropa/Kinematics.h"


namespace livpropa {

///////////////////////////////////////////////////////////////////////////////////////////////////

double AbstractKinematics::computeEnergyFromMomentum(const double& p, const int& id) const {
	return sqrt(computeEnergy2FromMomentum(p, id));
}

bool AbstractKinematics::isSpecialRelativistic() const {
	return getNameTag() == "SR";
}

bool AbstractKinematics::isLorentzViolatingMonochromatic() const {
	if (getNameTag().find("LIVmono") != std::string::npos)
		return true;

	return false;
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

}

double LorentzViolatingKinematics::selectFinalMomentumSmallest(const vector<double>& ps) {
	return *std::min_element(ps.begin(), ps.end());
}

double LorentzViolatingKinematics::selectFinalMomentumLargest(const vector<double>& ps) {
	return *std::max_element(ps.begin(), ps.end());
}

double LorentzViolatingKinematics::selectFinalMomentumAverage(const vector<double>& ps) {
	double pTot = 0;
	for (auto& p : ps)
		pTot += p;

	return pTot / ps.size();
}

double LorentzViolatingKinematics::selectFinalMomentumClosest(const vector<double>& ps, const double& E, const int& id) {
	SpecialRelativisticKinematics* sr = new SpecialRelativisticKinematics();
	double p0 = sr->computeMomentumFromEnergy(E, id);
	
	vector<double> dps;
	for (auto& p : ps)
		dps.push_back(p - p0);

	vector<double>::iterator idx = std::min_element(dps.begin(), dps.end());

	return dps[*idx];
}

string LorentzViolatingKinematics::info() const {
	return "LorentzViolatingKinematics";
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template<typename LIV>
LorentzViolatingKinematicsMonochromatic<LIV>::LorentzViolatingKinematicsMonochromatic(double chi, SymmetryBreaking symmetryBreaking) {
	setCoefficient(chi);
	setSymmetryBreaking(symmetryBreaking);
}

template<typename LIV>
LorentzViolatingKinematicsMonochromatic<LIV>::~LorentzViolatingKinematicsMonochromatic() {
}

template<typename LIV>
void LorentzViolatingKinematicsMonochromatic<LIV>::setOrder(int n) {
	order = n;
}

template<typename LIV>
void LorentzViolatingKinematicsMonochromatic<LIV>::setCoefficient(double coeff) {
	coefficient = coeff;
}

template<typename LIV>
int LorentzViolatingKinematicsMonochromatic<LIV>::getOrder() const {
	return static_cast<const LIV*>(this)->getOrder();
}

template<typename LIV>
double LorentzViolatingKinematicsMonochromatic<LIV>::getCoefficient() const {
	return coefficient;
}

template<typename LIV>
string LorentzViolatingKinematicsMonochromatic<LIV>::getNameTag() const {
	return "LIVmono" + std::to_string(order);
}

template<typename LIV>
string LorentzViolatingKinematicsMonochromatic<LIV>::getIdentifier() const {
	char s[64];
	sprintf(s, "LIV%i_chi_%+2.1e", order, coefficient);
	return string(s);
}

template<typename LIV>
double LorentzViolatingKinematicsMonochromatic<LIV>::getSymmetryBreakingShift(const double& p) const {
	return coefficient * pow(p * c_light / energy_planck, order) * pow_integer<2>(p * c_light);
}

template<typename LIV>
string LorentzViolatingKinematicsMonochromatic<LIV>::info() const {
	string s = "";
	s += "LorentzViolatingKinematicsMonochromatic";
	s += "(n = " + std::to_string(order) + "; chi = " + std::to_string(coefficient) + ")"; 
	return s;
}

template<typename LIV>
double LorentzViolatingKinematicsMonochromatic<LIV>::computeEnergy2FromMomentum(const double& p, const int& id) const {
	double m = particleMasses.at(id);
	double ds = this->getSymmetryBreakingShift(p);
	return pow_integer<2>(p * c_light) + pow_integer<2>(m * c_squared) + ds;
}

template<typename LIV>
double LorentzViolatingKinematicsMonochromatic<LIV>::computeMomentumFromEnergy(const double& E, const int& id) const {
	auto kin = static_cast<const LIV*>(this);
	return kin->computeMomentumFromEnergy(E, id);
}

template<int N>
int LIVKinematicsMonochromatic<N>::getOrder() const {
	return N;
}

template<int N>
LIVKinematicsMonochromatic<N>::LIVKinematicsMonochromatic(double chi, SymmetryBreaking symmetryBreaking) : LorentzViolatingKinematicsMonochromatic<LIVKinematicsMonochromatic<N> >(0, chi, symmetryBreaking) {
}

template<>
double LIVKinematicsMonochromatic<0>::_computeMomentumFromEnergy(const double& E, const int& id) const {
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
double LIVKinematicsMonochromatic<1>::_computeMomentumFromEnergy(const double& E, const int& id) const {
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
double LIVKinematicsMonochromatic<2>::_computeMomentumFromEnergy(const double& E, const int& id) const {
	double m = particleMasses.at(id);
	double chi = coefficient;
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

	return selectFinalMomentum(ps, E, id);
}

template<int N>
double LIVKinematicsMonochromatic<N>::_computeMomentumFromEnergy(const double& E, const int& id) const {
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


///////////////////////////////////////////////////////////////////////////////////////////////////

Kinematics::Kinematics() {
	 specialRelativity = new SpecialRelativisticKinematics();
}

Kinematics::Kinematics(vector<int> p, vector<ref_ptr<AbstractKinematics>> kin) {
	if (p.size() != kin.size())
		throw std::length_error("Vector of particles and kinematics must have the same length.");

	for (size_t i = 0; i < p.size(); i++) 
		add(p[i], kin[i]);

	specialRelativity = new SpecialRelativisticKinematics();
}

Kinematics::Kinematics(vector<int> p, ref_ptr<AbstractKinematics> kin) {
	for (size_t i = 0; i < p.size(); i++) 
		add(p[i], kin);

	specialRelativity = new SpecialRelativisticKinematics();
}

Kinematics::Kinematics(vector<pair<int, ref_ptr<AbstractKinematics>>> kin) {
	for (size_t i = 0; i < kin.size(); i++) 
		add(kin[i].first, kin[i].second);

	specialRelativity = new SpecialRelativisticKinematics();
}

void Kinematics::add(const int& particle, const ref_ptr<AbstractKinematics>& kin) {
	kinematics[particle] = kin;
}

void Kinematics::remove(const int& particle) {
	if (! exists(particle))
		throw runtime_error("Cannot retrieve inexistent particle with id " + std::to_string(particle) + ".");

	kinematics.erase(particle);
}

bool Kinematics::isLorentzInvariant() const {
	for (auto& kin : kinematics) {
		if (kin.second->getNameTag() != "SR")
			return false;
	}
	return true;
}

bool Kinematics::isLorentzViolatingKinematics() const {
	return ! isLorentzInvariant();
}

bool Kinematics::exists(const int& pId) const {
	vector<int> particles = getParticles();
	return std::find(particles.begin(), particles.end(), pId) != particles.end();
}

vector<int> Kinematics::getParticles() const {
	vector<int> particles;
	for (auto& kin : kinematics)
		particles.push_back(kin.first);

	return particles;
}

string Kinematics::getIdentifierForParticle(const int& pId, bool showParticleId) const {
	char identifier[128] = "";
	string info = find(pId, false)->getFilenamePart().c_str();

	size_t s = 0;
	if (showParticleId)
		s += sprintf(identifier + s, "Id_%+i-%s", pId, info.c_str());
	else
		s += sprintf(identifier + s, "%s", info.c_str());

	return std::string(identifier);
}

string Kinematics::getIdentifier(const std::vector<int>& particles, bool simplify) const {
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

ParticleKinematicsMap Kinematics::getParticleKinematicsMap() const {
	return kinematics;
}

const ref_ptr<AbstractKinematics>& Kinematics::find(const int& id, bool showWarningInexistent) const {
	bool declared = exists(id);

	if (declared) {
		ParticleKinematicsIterator it = kinematics.find(id);
		return (*it).second;
	} else {
		if (showWarningInexistent) 
			KISS_LOG_WARNING << "Cannot retrieve inexistent particle with id " << id << "." << "Returning special-relativistic kinematics." << endl;
		return specialRelativity;
	}		
}

const ref_ptr<AbstractKinematics>& Kinematics::operator[](const int& pId) {
	return find(pId);
}

ref_ptr<AbstractKinematics> Kinematics::operator[](const int& pId) const {
	return find(pId);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template<int N>
LorentzViolatingKinematicsMonochromatic<LIVKinematicsMonochromatic<N>> convertToLorentzViolatingKinematicsMonochromatic(const SpecialRelativisticKinematics& kin) {
	LorentzViolatingKinematicsMonochromatic<LIVKinematicsMonochromatic<N>>* k = new LorentzViolatingKinematicsMonochromatic<LIVKinematicsMonochromatic<N>>(0.);
	return *k;
}

std::ostream& operator<<(std::ostream& os, const AbstractKinematics& kin) {
	os << kin.info();
	return os;
}

std::ostream& operator<<(std::ostream& os, const Kinematics& kin) {
	os << "Kinematics: " << endl;

	vector<int> particles = kin.getParticles();
	for (auto& p : particles) 
		os << "  . particle " << p << " ==> " << kin[p]->info() << endl;

	return os;
}



// explicit instantiations for required template parameters
template class LorentzViolatingKinematicsMonochromatic<LIVKinematicsMonochromatic<0>>;
template class LorentzViolatingKinematicsMonochromatic<LIVKinematicsMonochromatic<1>>;
template class LorentzViolatingKinematicsMonochromatic<LIVKinematicsMonochromatic<2>>;


} // namespace livpropa