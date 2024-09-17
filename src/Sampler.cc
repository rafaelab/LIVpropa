#include "livpropa/Sampler.h"


namespace livpropa {


/****************************************************************************/

SamplerEventsEnergy::SamplerEventsEnergy(int pId, double s) {
	setSampling(s);
	setParticleId(pId);
}

void SamplerEventsEnergy::setSampling(double s) {
	sampling = s;
}

void SamplerEventsEnergy::setParticleId(int pId) {
	particleId = pId;
}

double SamplerEventsEnergy::getSampling() const {
	return sampling;
}

int SamplerEventsEnergy::getParticleId() const {
	return particleId;
}

double  SamplerEventsEnergy::computeWeight(int id, double E, double f, unsigned int counter) const {
	if (id != particleId)
		return 0;

	if (sampling >= 1) {
		return 1;
	} else if (sampling < 0) {
		return 0;
	} else {
		Random& random = Random::instance();
		double r = weightFunction(id, E, f, counter);
		if (random.rand() < r) // accept and return weight
			return 1. / r;
		else // reject
			return 0;
	}

}


/****************************************************************************/

SamplerEventsUniform::SamplerEventsUniform(int pId, double s) {
	setSampling(s);
	setParticleId(pId);
}

void SamplerEventsUniform::setSampling(double s) {
	sampling = s;
}

void SamplerEventsUniform::setParticleId(int pId) {
	particleId = pId;
}

double SamplerEventsUniform::getSampling() const {
	return sampling;
}

int SamplerEventsUniform::getParticleId() const {
	return particleId;
}

double SamplerEventsUniform::computeWeight(int id, double E, double f, unsigned int counter) const {
	if (id != particleId)
		return 0;

	if (sampling >= 1.) {
		return 1.;
	} else if (sampling <= 0) {
		return 0.;
	} else {
		// acceptance-rejection method
		Random& random = Random::instance();
		return random.rand() < sampling ? 1. / sampling : 0.;
	}
}


/****************************************************************************/

SamplerEventsEnergyFractionPowerLaw::SamplerEventsEnergyFractionPowerLaw(double idx, int pId, double s) : SamplerEventsEnergy(pId, s) {
	setIndex(idx);
}

void SamplerEventsEnergyFractionPowerLaw::setIndex(double idx) {
	index = idx;
}

double SamplerEventsEnergyFractionPowerLaw::getIndex() const {
	return index;
}

double SamplerEventsEnergyFractionPowerLaw::weightFunction(int id, double E, double f, unsigned int counter) const {
	return pow(f, (1 - sampling) * index);
}



// /****************************************************************************/

// SamplerEventsEnergyNormal::SamplerEventsEnergyNormal(int pId, double s, double mu, double sigma) : SamplerEventsEnergy(pId, s) {
// 	setMean(mu);
// 	setStandardDeviation(sigma);
// }

// void SamplerEventsEnergyNormal::setMean(double mu) {
// 	mean = mu;
// }

// void SamplerEventsEnergyNormal::setStandardDeviation(double sigma) {
// 	standardDeviation = sigma;
// }

// double SamplerEventsEnergyNormal::getMean() const {
// 	return mean;
// }

// double SamplerEventsEnergyNormal::getStandardDeviation() const {
// 	return standardDeviation;
// }

// double SamplerEventsEnergyNormal::weightFunction(int id, double E, double f, int counter) const {
// 	return pow(exp(- pow_integer<2>(E - mean) / 2. / pow_integer<2>(standardDeviation)), 1 - sampling);
// }



// /****************************************************************************/

// SamplerEventsEnergyLogNormal::SamplerEventsEnergyLogNormal(int pId, double s, double mu, double sigma) : SamplerEventsEnergy(pId, s) {
// 	setMean(mu);
// 	setStandardDeviation(sigma);
// }

// void SamplerEventsEnergyLogNormal::setMean(double mu) {
// 	mean = mu;
// }

// void SamplerEventsEnergyLogNormal::setStandardDeviation(double sigma) {
// 	standardDeviation = sigma;
// }

// double SamplerEventsEnergyLogNormal::getMean() const {
// 	return mean;
// }

// double SamplerEventsEnergyLogNormal::getStandardDeviation() const {
// 	return standardDeviation;
// }

// double SamplerEventsEnergyLogNormal::weightFunction(int id, double E, double f, int counter) const {
// 	return pow(exp(- pow_integer<2>(log(E) - mean) / 2. / pow_integer<2>(standardDeviation)) / (E / mean), 1 - sampling);
// }



/*
***************************************************************************/
SamplerEventsNull::SamplerEventsNull() {
}

double SamplerEventsNull::computeWeight(int id, double E, double f, unsigned int counter) const {
	return 1.;
}

/****************************************************************************/

SamplerEventsList::SamplerEventsList() {
}

SamplerEventsList::SamplerEventsList(std::vector<ref_ptr<SamplerEvents>> s) {
	for (size_t i = 0; i < s.size(); i++) {
		add(s[i]);
	}
}

void SamplerEventsList::add(SamplerEvents *SamplerEvents) {
	samplers.push_back(SamplerEvents);
}

double SamplerEventsList::computeWeight(int id, double E, double f, unsigned int counter) const {
	if (samplers.size() == 0)
		return 1.;

	double w = 0;
	int k = 0;
	for (size_t i = 0; i < samplers.size(); i++) {
		double w0 = samplers[i]->computeWeight(id, E, f, counter);
		if (w0 > 0) {
			if (k == 0)
				w = 1;
			w *= w0;
			k++;
		}
	}

	return w;
}


/****************************************************************************/

SamplerDistributionUniform::SamplerDistributionUniform(double vmin, double vmax, unsigned int nBins, std::string scale) {
	ref_ptr<Histogram1D> h = new Histogram1D(vmin, vmax, nBins, scale);	
	setSize(0);
	setDistribution(h);
}

void SamplerDistributionUniform::initDistribution(double vmin, double vmax, unsigned int nBins, std::string scale) {
	ref_ptr<Histogram1D> h = new Histogram1D(vmin, vmax, nBins, scale);	
	setDistribution(h);
}

void SamplerDistributionUniform::setDistribution(ref_ptr<Histogram1D> dist) {
	distribution = dist;
}

void SamplerDistributionUniform::setSize(int n) {
	datasetSize = n;
}

ref_ptr<Histogram1D> SamplerDistributionUniform::getDistribution() const {
	return distribution;
}

unsigned int SamplerDistributionUniform::getSize() const {
	return datasetSize;
}

double SamplerDistributionUniform::getSample() const {
	return distribution->getSample();
}

double SamplerDistributionUniform::interpolateAt(const double &v) const {
	return distribution->interpolateAt(v);
}

void SamplerDistributionUniform::prepareCDF() {
	distribution->prepareCDF();
}

void SamplerDistributionUniform::append(const std::vector<double>& v) {
	for (size_t i = 0; i < v.size(); i++) {
		distribution->push(v[i]);
	}
	datasetSize = v.size();
}

void SamplerDistributionUniform::push(const double& v) {
	distribution->push(v);
	datasetSize++;
}

void SamplerDistributionUniform::clear() {
	distribution->clear();
	datasetSize = 0;
}


// /****************************************************************************/

// SamplerDistributionList::SamplerDistributionList() {
// }

// SamplerDistributionList::SamplerDistributionList(std::vector<ref_ptr<SamplerDistribution>> s) {
// 	for (size_t i = 0; i < s.size(); i++) {
// 		add(s[i]);
// 	}
// }

// void SamplerDistributionList::add(SamplerDistribution *SamplerDistribution) {
// 	samplers.push_back(SamplerEvents);
// }

// void SamplerDistributionList::transformToCDF() {
// 	for (size_t i = 0; i < samplers.size(); i++) {
// 		samplers[i]->transformToCDF();
// 	}
// }

// void SamplerDistributionList::clear() {
// 	for (size_t i = 0; i < samplers.size(); i++) {
// 		samplers[i]->clear();
// 	}
// }

// std::vector<double> SamplerDistributionList::getSample(int nSamples) const {
// 	// std::vector<std::vector<double>> sample;

// 	// for (size_t j = 0; j < samplers.size(); j++) {
// 	// 	std::vector<double> tmp;
// 	// 	for (size_t i = 0; i < std::min(nSamples, distribution->getNumberOfBins()); i++) {
// 	// 		tmp.push_
// 	// 	}
// 	// 	sample.push_back(distribution->getSample());
// 	// }

// 	// return sample;
// }


// /****************************************************************************/
// SamplerDistributionNull::SamplerDistributionNull() {
// }

// std::vector<double> SamplerDistributionNull::getSample(int nSamples) const {
// 	std::vector<double> v;
// 	for (size_t i = 0; i < nSamples; i++) {
// 		v.push_back(0.);
// 	}

// 	return v;
// }


// SamplerEventsBinnedDistribution::SamplerEventsBinnedDistribution(int pId, int nSamples, double vmin, double vmax, int nBins, std::string scale, double normalisation) {
// 	ref_ptr<Histogram1D<double, double>> h = new Histogram1D<double, double>(vmin, vmax, nBins, scale, normalisation);
// 	setDistribution(h);
// 	setParticleId(pId);
// 	setNumberOfSamples(nSamples);
// 	setGlobalCounter(0);
// }

// void SamplerEventsBinnedDistribution::setParticleId(int pId) {
// 	particleId = pId;
// }

// void SamplerEventsBinnedDistribution::setGlobalCounter(int i) {
// 	globalCounter = i;
// }

// void SamplerEventsBinnedDistribution::setStatus(bool b) {
// 	isHistogramFilled = b;
// }

// void SamplerEventsBinnedDistribution::setNumberOfSamples(int n) {
// 	nSamples = n;
// }

// void SamplerEventsBinnedDistribution::setDistribution(ref_ptr<Histogram1D<double, double>> h) {
// 	distribution = h;
// }

// bool SamplerEventsBinnedDistribution::getStatus() const {
// 	return isHistogramFilled;
// }

// ref_ptr<Histogram1D<double, double>> SamplerEventsBinnedDistribution::getDistribution() const {
// 	return distribution;
// }

// double SamplerEventsBinnedDistribution::computeWeight(int id, double E, double f, int counter) const {
// 	if (id != particleId)
// 		return 0.;

// 	if (globalCounter > nSamples)
// 		return 0.;
// 	globalCounter++;

// 	if (getStatus()) {
// 		// Random& random = Random::instance();
// 		// size_t bin = random.randBin(binContents);
// 		return 1;

// 	} else {
// 		return 1;
// 	}
// }

// void SamplerEventsBinnedDistribution::clear() {
// 	setGlobalCounter(0);
// 	distribution->clear();
// }



/****************************************************************************/

} // namespace crpropa