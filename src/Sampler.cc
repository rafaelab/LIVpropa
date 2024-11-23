#include "livpropa/Sampler.h"

namespace livpropa {


void Sampler::setHistogram(ref_ptr<AbstractHistogram1D> h) {
	histogram = h;
}

ref_ptr<AbstractHistogram1D> Sampler::getHistogram() const {
	return histogram;
}

void Sampler::computeCDF() {
	vector<double> edges = histogram->getBinEdges();
	vector<double> contents = histogram->getBinContents();
	size_t nBins = histogram->getNumberOfBins();

		if (edges.size() < 2)
		throw std::runtime_error("Number of bin edges must be greater than one.");

	// ensure vector is empty
	cdf.clear();

	// first entry of CDF is 0
	cdf.push_back(0);

	double cum = 0;
	for (size_t i = 1; i < nBins; i++) {
		auto bin = histogram->getBin(i);
		double dx = bin->getWidth();
		double y = contents[i];
		cum += y * dx;
		cdf.push_back(cum);
	}

	for (size_t i = 1; i < nBins; i++) {
		cdf[i] /= cdf.back();
	}
}



SamplerInverse::SamplerInverse() {
}

SamplerInverse::SamplerInverse(ref_ptr<AbstractHistogram1D> h) {
	setHistogram(h);
	computeCDF();
}

// SamplerInverse::~SamplerInverse() {
// }

double SamplerInverse::getSample(Random& random) const {
	double r = random.rand();

	auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
	size_t idx = std::distance(cdf.begin(), it) - 1;

	// if (idx >= edges.size() - 1)
	// 	idx = edges.size() - 2;

	return interpolate(r, cdf, histogram->getBinCentres());
}

vector<double> SamplerInverse::getSamples(Random& random, unsigned int nSamples) const {
	vector<double> samples;

	for (size_t i = 0; i < nSamples; i++) {
		samples.push_back(getSample(random));
	}

	return samples;
}


} // namespace livpropa
