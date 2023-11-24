#include "livpropa/Kinematics.h"



namespace livpropa {


double Kinematics::computeEnergy2FromMomentum(const double& p, const double& m) const {
	double ds = getSymmetryBreakingShift(p);
	return pow_integer<2>(p * c_light) + pow_integer<2>(m * c_squared) + ds;
}

double Kinematics::computeEnergyFromMomentum(const double& p, const double& m) const {
	return sqrt(computeEnergy2FromMomentum(p, m));
}



SpecialRelativity::SpecialRelativity() {
}

std::string SpecialRelativity::getShortIdentifier() const {
	return "SR";
}

std::string SpecialRelativity::getLocationData(std::vector<int> particles) const {
	return "";
}

double SpecialRelativity::getSymmetryBreakingShift(const double& p) const {
	return 0;
}

double SpecialRelativity::computeEnergyFromMomentum(const double& p, const double&m ) const {
	return hypot(p * c_light, m * c_squared);
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

double MonochromaticLIV::getSymmetryBreakingShift(const double& p) const {
	return 0;
}

double MonochromaticLIV::computeEnergyFromMomentum(const double& p, const double& m) const {
	return 0;
}



} // namespace livpropa