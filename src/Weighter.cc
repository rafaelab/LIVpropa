#include "livpropa/Weighter.h"



namespace livpropa {


void Weighter::setWeighterType(WeighterType t) {
	type = t;
}

WeighterType Weighter::getType() const {
	return type;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

WeighterUniformFraction::WeighterUniformFraction(int pId, double s) {
	setSamplingFraction(s);
	setParticleId(pId);
	setWeighterType(WeighterType::UniformFraction);
}

void WeighterUniformFraction::setSamplingFraction(double s) {
	samplingFraction = s;
}

void WeighterUniformFraction::setParticleId(int pId) {
	particleId = pId;
}

double WeighterUniformFraction::getSamplingFraction() const {
	return samplingFraction;
}

int WeighterUniformFraction::getParticleId() const {
	return particleId;
}

double WeighterUniformFraction::computeWeight(const int& id, const double& E, const double& f, const int& counter, Random& random) const {
	if (id != particleId)
		return 0;

	if (samplingFraction >= 1.) {
		return 1;
	} else if (samplingFraction <= 0) {
		return 0;
	} else {
		if (random.rand() < samplingFraction) // accept and return weight
			return 1. / samplingFraction;
		else // reject
			return 0;
	}
}

string WeighterUniformFraction::getNameTag() const {
	return "uniformFraction";
}



///////////////////////////////////////////////////////////////////////////////////////////////////

WeighterUniformNumber::WeighterUniformNumber(int pId, unsigned int n) {
	setNumberOfEvents(n);
	setParticleId(pId);
	setWeighterType(WeighterType::UniformNumber);
}

void WeighterUniformNumber::setNumberOfEvents(unsigned int n) {
	nEvents = n;
}

void WeighterUniformNumber::setParticleId(int pId) {
	particleId = pId;
}

unsigned int WeighterUniformNumber::getNumberOfEvents() const {
	return nEvents;
}

int WeighterUniformNumber::getParticleId() const {
	return particleId;
}

double WeighterUniformNumber::computeWeight(const int& id, const double& E, const double& f, const int& counter, Random& random) const {
	if (id != particleId)
		return 0;

	// UNTESTED!!!!!!!!!!!!

	if (counter >= nEvents)
		return 0;
	
	return 1;
}

string WeighterUniformNumber::getNameTag() const {
	return "uniformNumber";
}



///////////////////////////////////////////////////////////////////////////////////////////////////

WeighterList::WeighterList() {
	setWeighterType(WeighterType::List);
}

WeighterList::WeighterList(std::vector<ref_ptr<Weighter>> ws) {
	for (size_t i = 0; i < ws.size(); i++) {
		add(ws[i]);
	}
	setWeighterType(WeighterType::List);
}

void WeighterList::add(Weighter* w) {
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

string WeighterList::getNameTag() const {
	return "list";
}



///////////////////////////////////////////////////////////////////////////////////////////////////

WeighterNull::WeighterNull() {
	setWeighterType(WeighterType::Null);
}

double WeighterNull::computeWeight(const int& id, const double& E, const double& f, const int& counter, Random& random) const {
	return 1;
}

string WeighterNull::getNameTag() const {
	return "null";
}


} // namespace livpropa
