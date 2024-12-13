#include "livpropa/Sampler.h"

namespace livpropa {


string getSamplerNameTag(SamplerType t) {
	switch (t) {
		case SamplerType::Inverse:
			return "inverse";
		case SamplerType::Rejection:
			return "rejection";
		case SamplerType::Importance:
			return "importance";
		case SamplerType::Nested:
			return "nested";
		case SamplerType::MCMC:
			return "mcmc";
		default:
			throw std::runtime_error("Unknown sampler type.");
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void Sampler::setDistribution(ref_ptr<Histogram1D> h) {
	histogram = h;
}

ref_ptr<Histogram1D> Sampler::getDistribution() const {
	return histogram;
}

vector<std::pair<double, double>> Sampler::getSamples(unsigned int nSamples, Random& random, const std::pair<double, double>& range) const {
	vector<std::pair<double, double>> samples;
	for (unsigned int i = 0; i < nSamples; i++) {
		samples.push_back(getSample(random, range));
	}

	return samples;
}

void Sampler::setType(SamplerType t) {
	type = t;
}

SamplerType Sampler::getType() const {
	return type;
}

string Sampler::getNameTag() const {
	return getSamplerNameTag(type);
}

void Sampler::push(const double& v, const double& w) {
	histogram->push(v, w);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

InverseSampler::InverseSampler() {
	setType(SamplerType::Inverse);
}

InverseSampler::InverseSampler(ref_ptr<Histogram1D> h) {
	setType(SamplerType::Inverse);
	if (h->getIsPDF())
		setDistribution(h);
	else
		h->makePDF();
	setDistribution(h);
	computeCDF();
}

void InverseSampler::computeCDF() {
	cdf = histogram->computeVectorCDF();
}

std::pair<double, double> InverseSampler::getSample(Random& random, const std::pair<double, double>& range) const {
	double r = random.randUniform(range.first, range.second);
	unsigned int nBins = histogram->getNumberOfBins();

	if (r > cdf[nBins - 1])
		return std::make_pair(histogram->getBinCentre(0), 1.);

	double x = interpolate(r, cdf, histogram->getBinCentres());
	double w = 1.;

	return std::make_pair(x, w);
}

void InverseSampler::reset() {
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

void RejectionSampler::computeCDF() {
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

ref_ptr<Histogram1D> RejectionSampler::getCDF() const {
	ref_ptr<Histogram1D> hCum = histogram;
	hCum->setBinContents(cdf);
	return hCum;
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

void RejectionSampler::reset() {
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

	// apply consistent naming
	std::unordered_map<string, string> names = {
		{"constant", "uniform"},
		{"log", "log10"},
		{"exponential", "exp"},
		{"sqrt", "power0.5"},
		{"linear", "power1"},
		{"quadratic", "power2"},
		{"cubic", "power3"},
		{"inverselinear", "inversepower1"},
		{"inversequadratic", "inversepower2"},
		{"inversecubic", "inversepower3"},
		{"inversesqrt", "inversepower0.5"},
		{"inverseexp", "inverseexponential"}
	};
	if (names.find(weight) != names.end()) 
		weight = names[weight];


	if (weight == "uniform") {
		setWeightFunction([](const double& x) { return 1; });
	} else if (weight == "log10" ) {
		setWeightFunction([](const double& x) { return abs(log10(x)); });
	} else if (weight == "ln") {
		setWeightFunction([](const double& x) { return abs(log(x)); });
	} else if (weight == "exponential") {
		setWeightFunction([](const double& x) { return exp(x); });
	} else if (weight == "inverseexponential") {
		setWeightFunction([](const double& x) { return exp(- x); });
	} else if (weight.find("power") == 0) {
		double power = parseWeightFunctionName(weight, "power");
		setWeightFunction([power](const double& x) { return pow(x, power) / (1 + power); });
	} else if (weight.find("inversepower") == 0) {
		double power = parseWeightFunctionName(weight, "inversepower");
		setWeightFunction([power](const double& x) { return pow(x, - power) / (1 - power); });
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

void ImportanceSampler::computeCDF() {
	inverseSampler.computeCDF();
}

ref_ptr<Histogram1D> ImportanceSampler::getCDF() const {
	ref_ptr<Histogram1D> hCum = histogram;
	hCum->setBinContents(cdf);
	return hCum;
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

double ImportanceSampler::parseWeightFunctionName(const string& str, const string& pattern) {
	if (str.find("power") == 0) {
		string powerStr = str.substr(pattern.length());
		std::istringstream iss(powerStr);
		double power;
		iss >> power;
		return power;
	}
	return 1.;
}

void ImportanceSampler::reset() {
}



///////////////////////////////////////////////////////////////////////////////////////////////////

NestedSampler::NestedSampler(unsigned int nLivePoints, unsigned int maxIterations) {
	setType(SamplerType::Nested);
	setMaximumIterations(maxIterations);
	setNumberOfLivePoints(nLivePoints);
	logEvidence = -std::numeric_limits<double>::infinity();
	logWeight = 0;
}

NestedSampler::NestedSampler(ref_ptr<Histogram1D> h, unsigned int nLivePoints) {
	setType(SamplerType::Nested);
	setDistribution(h);
	setNumberOfLivePoints(nLivePoints);
	setLikelihoodFunction(h->getInterpolator());
	logEvidence = -std::numeric_limits<double>::infinity();
	logWeight = 0;
}

void NestedSampler::setMaximumIterations(unsigned int n) {
	maxIterations = n;
}

void NestedSampler::setNumberOfLivePoints(unsigned int n) {
	nLivePoints = n;
}

void NestedSampler::setLikelihoodFunction(std::function<double(double)> func) {
	likelihood = func;
}

void NestedSampler::setDistribution(ref_ptr<Histogram1D> h) {
	histogram = h;
	likelihood = h->getInterpolator();
}

unsigned int NestedSampler::getNumberOfLivePoints() const {
	return nLivePoints;
}

unsigned int NestedSampler::getMaximumIterations() const {
	return maxIterations;
}

std::function<double(double)> NestedSampler::getLikelihoodFunction() const {
	return likelihood;
}

std::pair<double, double> NestedSampler::getSample(Random& random, const std::pair<double, double>& range) const {
	// initialise live points (if not already done)
	if (livePoints.empty()) {
		for (unsigned int i = 0; i < nLivePoints; ++i) {
			double point = random.randUniform(range.first, range.second);
			livePoints.push_back(point);
			liveLikelihoods.push_back(likelihood(point));
		}
	}


	double newPoint = 0;
	double newLikelihood = 0;
	
	unsigned int iteration = 0;
	while (iteration < maxIterations) {
		// find the live point with the lowest likelihood
		auto minIt = std::min_element(liveLikelihoods.begin(), liveLikelihoods.end());
		size_t minIdx = std::distance(liveLikelihoods.begin(), minIt);
		double minLikelihood = *minIt;

		// update log evidence
		double logVolume = - static_cast<double>(minIdx + 1) / nLivePoints;
		logWeight = logVolume + minLikelihood;
		logEvidence = logSumExp(logEvidence, logWeight);

		// replace the worst live point with a new sample
		do {
			newPoint = random.randUniform(range.first, range.second);
			newLikelihood = likelihood(newPoint);
		} while (newLikelihood <= minLikelihood);

		livePoints[minIdx] = newPoint;
		liveLikelihoods[minIdx] = newLikelihood;

		iteration++;
	}

	return std::make_pair(newPoint, newLikelihood);
}

double NestedSampler::getLogEvidence() const {
	return logEvidence;
}

double NestedSampler::logSumExp(double a, double b) const {
	if (a == - std::numeric_limits<double>::infinity()) 
		return b;
	
	if (b == - std::numeric_limits<double>::infinity()) 
		return a;

	if (a > b) 
		return a + log1p(exp(b - a));
	
	return b + log1p(exp(a - b));
}

void NestedSampler::reset() {
	livePoints.clear();
	liveLikelihoods.clear();
	logEvidence = - std::numeric_limits<double>::infinity();
	logWeight = 0;
}



} // namespace livpropa
