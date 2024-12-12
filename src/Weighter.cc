#include "livpropa/Weighter.h"



namespace livpropa {


void RuntimeWeighter::setType(Type t) {
	type = t;
}

RuntimeWeighter::Type RuntimeWeighter::getType() const {
	return type;
}

string RuntimeWeighter::getNameTag() const {
	switch (type) {
		case Type::Null:
			return "null";
		case Type::List:
			return "list";
		case Type::EnergyFraction:
			return "energy fraction";
		case Type::EnergyFractionUniform:
			return "energy fraction: uniform";
		case Type::EnergyFractionPowerLaw:
			return "energy fraction: power law";
		default:
			return "unknown";
	}
}

void RuntimeWeighter::reset() {
}


///////////////////////////////////////////////////////////////////////////////////////////////////

WeighterEnergyFraction::WeighterEnergyFraction() {
	setType(RuntimeWeighter::Type::EnergyFraction);
}

WeighterEnergyFraction::WeighterEnergyFraction(int pId) {
	setType(RuntimeWeighter::Type::EnergyFraction);
	setParticleId(pId);
}

WeighterEnergyFraction::WeighterEnergyFraction(int pId, std::function<double(double)> f) {
	setType(RuntimeWeighter::Type::EnergyFraction);
	setParticleId(pId);
	setWeightFunction(f);
}

WeighterEnergyFraction::~WeighterEnergyFraction() {
}

void WeighterEnergyFraction::setParticleId(int pId) {
	particleId = pId;
}

void WeighterEnergyFraction::setWeightFunction(std::function<double(double)> f) {
	weightFunction = f;
}

int WeighterEnergyFraction::getParticleId() const {
	return particleId;
}

std::function<double(double)> WeighterEnergyFraction::getWeightFunction() const {
	return weightFunction;
}

double WeighterEnergyFraction::computeWeight(const int& id, const double& E, const double& f, const int& counter, Random& random) const {
	if (id != particleId)
		return 0;
	
	// accept/reject and return weight
	double w = weightFunction(f);
	return (random.rand() < w) ? 1. / w : 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////

WeighterEnergyFractionUniform::WeighterEnergyFractionUniform() : WeighterEnergyFraction() {
	setType(RuntimeWeighter::Type::EnergyFractionUniform);
}

WeighterEnergyFractionUniform::WeighterEnergyFractionUniform(int pId, double s) : WeighterEnergyFraction(pId) {
	setType(RuntimeWeighter::Type::EnergyFractionUniform);
	setSamplingFraction(s);
	setParticleId(pId);

	std::function<double(double)> f = [s](const double& x) { 
		return s; 
		};
	setWeightFunction(f);
}

void WeighterEnergyFractionUniform::setSamplingFraction(double s) {
	samplingFraction = s;
}

double WeighterEnergyFractionUniform::getSamplingFraction() const {
	return samplingFraction;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

WeighterEnergyFractionPowerLaw::WeighterEnergyFractionPowerLaw() : WeighterEnergyFraction() {
	setType(RuntimeWeighter::Type::EnergyFractionPowerLaw);
}

WeighterEnergyFractionPowerLaw::WeighterEnergyFractionPowerLaw(int pId, double s) : WeighterEnergyFraction(pId) {
	setType(RuntimeWeighter::Type::EnergyFractionPowerLaw);
	setExponent(s);
	setParticleId(pId);

	std::function<double(double)> f = [s](const double& x) { 
		return pow(x, s); 
		};
	setWeightFunction(f);
}

void WeighterEnergyFractionPowerLaw::setExponent(double s) {
	exponent = s;
}

double WeighterEnergyFractionPowerLaw::getExponent() const {
	return exponent;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

WeighterList::WeighterList() {
	setType(Type::List);
}

WeighterList::WeighterList(std::vector<ref_ptr<RuntimeWeighter>> ws) {
	for (size_t i = 0; i < ws.size(); i++) {
		add(ws[i]);
	}
	setType(Type::List);
}

void WeighterList::add(RuntimeWeighter* w) {
	weighters.push_back(w);
}

double WeighterList::computeWeight(const int& id, const double& E, const double& f, const int& counter, Random& random) const {
	if (weighters.size() == 0)
		return 1.;

	double w = 0;
	int k = 0;
	for (size_t i = 0; i < weighters.size(); i++) {
		double w0 = weighters[i]->computeWeight(id, E, f, counter, random);
		if (w0 > 0) {
			if (k == 0)
				w = 1;
			w *= w0;
			k++;
		}
	}

	return w;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

WeighterNull::WeighterNull() {
	setType(Type::Null);
}

double WeighterNull::computeWeight(const int& id, const double& E, const double& f, const int& counter, Random& random) const {
	return 1;
}



///////////////////////////////////////////////////////////////////////////////////////////////////

void PosterioriWeighter::setType(Type t) {
	type = t;
}

PosterioriWeighter::Type PosterioriWeighter::getType() const {
	return type;
}

string PosterioriWeighter::getNameTag() const {
	switch (type) {
		case Type::Distribution:
			return "distribution";
		default:
			return "unknown";
	}
}

void PosterioriWeighter::reset() {
}


///////////////////////////////////////////////////////////////////////////////////////////////////

WeighterDistribution::WeighterDistribution() {
	setType(Type::Distribution);
}

WeighterDistribution::WeighterDistribution(int pId, unsigned int n) {
	setType(Type::Distribution);
	setNumberOfEvents(n);
	setParticleId(pId);
}

WeighterDistribution::WeighterDistribution(int pId, ref_ptr<Histogram1D> h, unsigned int n) {
	setType(Type::Distribution);
	setHistogram(h);
	setNumberOfEvents(n);
	setParticleId(pId);
}

void WeighterDistribution::setNumberOfEvents(unsigned int n) {
	nEvents = n;
}

void WeighterDistribution::setParticleId(int pId) {
	particleId = pId;
}

void WeighterDistribution::setHistogram(ref_ptr<Histogram1D> h) {
	histogram = h;
}

void WeighterDistribution::push(const double& v, const double& w) {
	histogram->push(v, w);
}

unsigned int WeighterDistribution::getNumberOfEvents() const {
	return nEvents;
}

int WeighterDistribution::getParticleId() const {
	return particleId;
}

ref_ptr<Histogram1D> WeighterDistribution::getHistogram() const {
	return histogram;
}

double WeighterDistribution::computeWeight(const int& id, const double& E, const double& f, const int& counter, Random& random) const {
	if (id != particleId)
		return 0;

	if (counter >= nEvents)
		return 0;

	histogram->push(E);
	nEntries++;
	
	return 1;
}

void WeighterDistribution::reset() {
	histogram->reset();
	nEntries = 0;
}





} // namespace livpropa
