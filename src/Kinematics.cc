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


MonochromaticLIV::MonochromaticLIV() {
}

MonochromaticLIV::MonochromaticLIV(unsigned int n) {
	setOrder(n);
}

MonochromaticLIV::MonochromaticLIV(unsigned int order, std::unordered_map<int, double> coeffs) {
	setCoefficients(coeffs);
}

MonochromaticLIV::MonochromaticLIV(unsigned int order, std::vector<int> particles, std::vector<double> coeffs) {
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

	// select the largest among the positive roots
	// this is completely arbitrary and should be rethought
	std::vector<double> ps;
	for (auto& root : roots) {
		if (root.real() > 0.)
			ps.push_back(root.real());
	}
	std::sort(ps.begin(), ps.end(), std::greater<double>()); 
	
	return ps[0];
}



} // namespace livpropa