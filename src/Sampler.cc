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

void Sampler::setType(SamplerType t) {
	type = t;
}

SamplerType Sampler::getType() const {
	return type;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

InverseSampler::InverseSampler() {
	setType(SamplerType::Inverse);
}

InverseSampler::InverseSampler(ref_ptr<Histogram1D> h) {
	setDistribution(h);
	computeCDF();
	setType(SamplerType::Inverse);
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
	setType(SamplerType::Rejection);
}

RejectionSampler::RejectionSampler(ref_ptr<Histogram1D> h, ref_ptr<Histogram1D> proposal, double maxRatio) {
	setDistribution(h);
	setProposalPDF(proposal);
	computeCDF();
	setType(SamplerType::Rejection);
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
	setType(SamplerType::Importance);
}

ImportanceSampler::ImportanceSampler(std::function<double(double)> weight) {
	setWeightFunction(weight);
	setType(SamplerType::Importance);
}

ImportanceSampler::ImportanceSampler(string weight) {
	setType(SamplerType::Importance);

	auto parsePower = [](const std::string& str) -> double {
		if (str.find("power") == 0) {
			std::string powerStr = str.substr(5);
			std::istringstream iss(powerStr);
			double power;
			iss >> power;
			return power;
		}
		return 1.;
	};

	auto parseInversePower = [](const std::string& str) -> double {
		if (str.find("ipower") == 0) {
			std::string powerStr = str.substr(6);
			std::istringstream iss(powerStr);
			double power;
			iss >> power;
			return power;
		}
		return 1.;
	};

	if (weight == "constant" or weight == "uniform") {
		setWeightFunction([](const double& x) { return 1; });
	} else if (weight == "linear" or weight == "lin") {
		setWeightFunction([](const double& x) { return x; });
	} else if (weight == "quadratic" or weight == "quad") {
		setWeightFunction([](const double& x) { return x * x; });
	} else if (weight == "cubic") {
		setWeightFunction([](const double& x) { return x * x * x; });
	} else if (weight == "exponential" or weight == "exp") {
		setWeightFunction([](const double& x) { return exp(x); });
	} else if (weight == "sqrt") {
		setWeightFunction([](const double& x) { return sqrt(x); });
	} else if (weight.find("power") == 0) {
		double power = parsePower(weight);
		setWeightFunction([power](const double& x) { return pow(x, power); });
	} else if (weight == "ilinear" or weight == "ilin" or weight == "i1" or weight == "inverselinear") {
		setWeightFunction([](const double& x) { return 1. / x; });
	} else if (weight == "iquadratic" or weight == "iquad" or weight == "inversequadratic") {
		setWeightFunction([](const double& x) { return 1. / (x * x); });
	} else if (weight == "icubic" or weight == "inversecubic") {
		setWeightFunction([](const double& x) { return 1. / (x * x * x); });
	} else if (weight == "iexponential" or weight == "iexp" or weight == "inverseexponential" or weight == "inverseexp") {
		setWeightFunction([](const double& x) { return exp(-x); });
	} else if (weight == "isqrt" || weight == "inversesqrt") {
		setWeightFunction([](const double& x) { return 1. / sqrt(x); });
	} else if (weight.find("ipower") == 0) {
		double power = parseInversePower(weight);
		setWeightFunction([power](const double& x) { return 1. / pow(x, power); });
	} else {
		throw std::runtime_error("Unknown weight function. Try setting it manually.");
	}
}

ImportanceSampler::ImportanceSampler(ref_ptr<Histogram1D> h, ref_ptr<Histogram1D> proposal, std::function<double(double)> weight) {
	setDistribution(h);
	setProposalPDF(proposal);
	computeCDF();
	setWeightFunction(weight);
	setType(SamplerType::Importance);
}

ImportanceSampler::ImportanceSampler(ref_ptr<Sampler> sampler) {
	setType(SamplerType::Importance);
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
