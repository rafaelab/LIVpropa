#include "livpropa/Sampler.h"

namespace livpropa {


///////////////////////////////////////////////////////////////////////////////////////////////////

void Sampler::setDistribution(ref_ptr<Histogram1D> h) {
	histogram = h;
}

ref_ptr<Histogram1D> Sampler::getDistribution() const {
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

vector<std::pair<double, double>> Sampler::getSamples(unsigned int nSamples, Random& random, const std::pair<double, double>& range) const {
	vector<std::pair<double, double>> samples;
	for (unsigned int i = 0; i < nSamples; i++) {
		samples.push_back(getSample(random, range));
	}

	return samples;
}

ref_ptr<Histogram1D> Sampler::getCumulativeDistribution() const {
	ref_ptr<Histogram1D> hCum = histogram;
	hCum->setBinContents(cdf);
	return hCum;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

InverseSampler::InverseSampler() {
}

InverseSampler::InverseSampler(ref_ptr<Histogram1D> h) {
	setDistribution(h);
	computeCDF();
}

std::pair<double, double> InverseSampler::getSample(Random& random, const std::pair<double, double>& range) const {
	double r = random.randUniform(range.first, range.second);
	unsigned int nBins = histogram->getNumberOfBins();

	auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
	size_t idx = std::distance(cdf.begin(), it) - 1;

	if (idx >= nBins - 1)
		idx = nBins - 2;

	double x = interpolate(r, cdf, histogram->getBinCentres());
	double w = 1.;

	return std::make_pair(x, w);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

ImportanceSampler::ImportanceSampler() {
}

ImportanceSampler::ImportanceSampler(ref_ptr<Histogram1D> h, ref_ptr<Histogram1D> proposal) {
	setDistribution(h);
	setProposalPDF(proposal);
	computeCDF();
}

void ImportanceSampler::computeCDF() {
	inverseSampler.computeCDF();
}

void ImportanceSampler::setProposalPDF(ref_ptr<Histogram1D> proposal) {
	proposalPDF = proposal;
	inverseSampler.setDistribution(proposal);
	inverseSampler.computeCDF();
}

ref_ptr<Histogram1D> ImportanceSampler::getProposalPDF() const {
	return proposalPDF;
}

std::pair<double, double> ImportanceSampler::getSample(Random& random, const std::pair<double, double>& range) const {
	std::pair<double, double> sample = inverseSampler.getSample(random, range);
	double x = sample.first;
	double w = histogram->getBinContent(histogram->getBinIndex(x)) / proposalPDF->getBinContent(proposalPDF->getBinIndex(x));
	
	return std::make_pair(x, w);
}





} // namespace livpropa
