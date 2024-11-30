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

string InverseSampler::getNameTag() const {
	return "inverse";
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

RejectionSampler::RejectionSampler() {
}

RejectionSampler::RejectionSampler(ref_ptr<Histogram1D> h, ref_ptr<Histogram1D> proposal, double maxRatio) {
	setDistribution(h);
	setProposalPDF(proposal);
	computeCDF();
	this->maxRatio = maxRatio;
}

string RejectionSampler::getNameTag() const {
	return "rejection";
}

void RejectionSampler::computeCDF() {
	inverseSampler.computeCDF();
}

void RejectionSampler::setProposalPDF(ref_ptr<Histogram1D> proposal) {
	proposalPDF = proposal;
}

ref_ptr<Histogram1D> RejectionSampler::getProposalPDF() const {
	return proposalPDF;
}

std::pair<double, double> RejectionSampler::getSample(Random& random, const std::pair<double, double>& range) const {
	while (true) {
		std::pair<double, double> sample = inverseSampler.getSample(random, range);
		double x = sample.first;
		double u = random.randUniform(0, 1);
		double pdfValue = histogram->getBinContent(histogram->getBinIndex(x));
		double proposalValue = proposalPDF->getBinContent(proposalPDF->getBinIndex(x));
		if (u < pdfValue / (proposalValue * maxRatio)) {
			double weight = pdfValue / proposalValue;
			return std::make_pair(x, weight);
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////

ImportanceSampler::ImportanceSampler() {
}

ImportanceSampler::ImportanceSampler(std::function<double(double)> weight) {
	setWeightFunction(weight);
}

ImportanceSampler::ImportanceSampler(string weight) {
	if (weight == "constant" or weight == "uniform") {
		setWeightFunction([](const double& x) { return 1; });
	} else if (weight == "linear" or weight == "lin" or weight == "1") {
		setWeightFunction([](const double& x) { return x; });
	} else if (weight == "quadratic" or weight == "quad" or weight == "2") {
		setWeightFunction([](const double& x) { return x * x; });
	} else if (weight == "cubic" or weight == "3") {
		setWeightFunction([](const double& x) { return x * x * x; });
	} else if (weight == "exponential" or weight == "exp") {
		setWeightFunction([](const double& x) { return exp(x); });
	} else if (weight == "sqrt") {
		setWeightFunction([](const double& x) { return sqrt(x); });
	} else if (weight == "power0.001") {
		setWeightFunction([](const double& x) { return pow(x, 0.001); });
	} else if (weight == "power0.01") {
		setWeightFunction([](const double& x) { return pow(x, 0.01); });
	} else if (weight == "power0.1") {
		setWeightFunction([](const double& x) { return pow(x, 0.1); });
	} else if (weight == "power0.2") {
		setWeightFunction([](const double& x) { return pow(x, 0.2); });
	} else if (weight == "power0.3") {
		setWeightFunction([](const double& x) { return pow(x, 0.3); });
	} else if (weight == "power0.4") {
		setWeightFunction([](const double& x) { return pow(x, 0.4); });
	} else if (weight == "power0.5") {
		setWeightFunction([](const double& x) { return pow(x, 0.5); });
	} else if (weight == "power0.6") {
		setWeightFunction([](const double& x) { return pow(x, 0.6); });
	} else if (weight == "power0.7") {
		setWeightFunction([](const double& x) { return pow(x, 0.7); });
	} else if (weight == "power0.8") {
		setWeightFunction([](const double& x) { return pow(x, 0.8); });
	} else if (weight == "power0.9") {
		setWeightFunction([](const double& x) { return pow(x, 0.9); });
	} else if (weight == "power0.99") {
		setWeightFunction([](const double& x) { return pow(x, 0.99); });
	} else if (weight == "power0.999") {
		setWeightFunction([](const double& x) { return pow(x, 0.999); });
	} else if (weight == "power1.001") {
		setWeightFunction([](const double& x) { return pow(x, 1.001); });
	} else if (weight == "power1.01") {
		setWeightFunction([](const double& x) { return pow(x, 1.01); });
	} else if (weight == "power1.1") {
		setWeightFunction([](const double& x) { return pow(x, 1.1); });
	} else if (weight == "power1.2") {
		setWeightFunction([](const double& x) { return pow(x, 1.2); });
	} else if (weight == "power1.3") {
		setWeightFunction([](const double& x) { return pow(x, 1.3); });
	} else if (weight == "power1.4") {
		setWeightFunction([](const double& x) { return pow(x, 1.4); });
	} else if (weight == "power1.5") {
		setWeightFunction([](const double& x) { return pow(x, 1.5); });
	} else if (weight == "ilinear" or weight == "ilin" or weight == "i1" or weight == "inverselinear" or weight == "inverse1") {
		setWeightFunction([](const double& x) { return 1. / x; });
	} else if (weight == "iquadratic" or weight == "iquad" or weight == "i2" or weight == "inversequadratic" or weight == "inverse2") {
		setWeightFunction([](const double& x) { return 1. / (x * x); });
	} else if (weight == "icubic" or weight == "i3" or weight == "inversecubic" or weight == "inverse3") {
		setWeightFunction([](const double& x) { return 1. / (x * x * x); });
	} else if (weight == "iexponential" or weight == "iexp" or weight == "inverseexponential" or weight == "inverseexp") {
		setWeightFunction([](const double& x) { return exp(- x); });
	} else if (weight == "isqrt" or weight == "isqrt" or weight == "inversesqrt") {
		setWeightFunction([](const double& x) { return 1. / sqrt(x); });
	} else if (weight == "ipower0001" or weight == "ipower0.001" or weight == "inversepower0.001") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.001); });
	} else if (weight == "ipower0.01" or weight == "inversepower0.01") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.01); });
	} else if (weight == "ipower0.1" or  weight == "inversepower0.1") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.1); });
	} else if (weight == "ipower0.2" or  weight == "inversepower0.2") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.2); });
	} else if (weight == "ipower0.3" or weight == "inversepower0.3") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.3); });
	} else if (weight == "ipower0.4" or weight == "inversepower0.4") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.4); });
	} else if (weight == "ipower0.5" or weight == "inversepower0.5") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.5); });
	} else if (weight == "ipower0.6" or weight == "inversepower0.6") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.6); });
	} else if (weight == "ipower0.7" or weight == "inversepower0.7") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.7); });
	} else if (weight == "ipower0.8" or weight == "inversepower0.8") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.8); });
	} else if (weight == "ipower0.9" or weight == "inversepower0.9") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.9); });
	} else if (weight == "ipower0.99" or weight == "inversepower0.99") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.99); });
	} else if (weight == "ipower0.999" or weight == "inversepower0.999") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 0.999); });
	} else if (weight == "ipower0.999" or weight == "inversepower1.001") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 1.001); });
	} else if (weight == "ipower0.999" or weight == "inversepower1.01") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 1.01); });
	} else if (weight == "ipower1.1" or weight == "inversepower1.1") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 1.1); });
	} else if (weight == "ipower1.2" or weight == "inversepower1.2") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 1.2); });
	} else if (weight == "ipower1.3" or weight == "inversepower1.3") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 1.3); });
	} else if (weight == "ipower1.4" or weight == "inversepower1.4") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 1.4); });
	} else if (weight == "ipower1.5" or weight == "inversepower1.5") {
		setWeightFunction([](const double& x) { return 1. / pow(x, 1.5); });
	} else {
		throw std::runtime_error("Unknown weight function. Try setting it manually.");
	}
}

ImportanceSampler::ImportanceSampler(ref_ptr<Histogram1D> h, ref_ptr<Histogram1D> proposal, std::function<double(double)> weight) {
	setDistribution(h);
	setProposalPDF(proposal);
	computeCDF();
	setWeightFunction(weight);
}

ImportanceSampler::ImportanceSampler(ref_ptr<Sampler> sampler) {
	if (sampler->getNameTag() == "importance") {
		auto s = static_cast<ImportanceSampler*>(sampler.get());
		setDistribution(s->getDistribution());
		setProposalPDF(s->getProposalPDF());
		setWeightFunction(s->getWeightFunction());
		computeCDF();
	} else {
		throw("Cannot cast ref_ptr<Sampler> to ImportanceSampler.");
	}
}

string ImportanceSampler::getNameTag() const {
	return "importance";
}

void ImportanceSampler::computeCDF() {
	inverseSampler.computeCDF();
}

void ImportanceSampler::setProposalPDF(ref_ptr<Histogram1D> proposal) {
	proposalPDF = proposal;
	inverseSampler.setDistribution(proposal);
	inverseSampler.computeCDF();
}

void ImportanceSampler::setWeightFunction(std::function<double(double)> func) {
	weightFunction = func;
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

std::function<double(double)> ImportanceSampler::getWeightFunction() const {
	return weightFunction;
}

} // namespace livpropa
