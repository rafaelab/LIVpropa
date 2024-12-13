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
		setWeightFunction([power](const double& x) { return 1. / pow(x, power); });
	} else if (weight.find("inversepower") == 0) {
		double power = parseWeightFunctionName(weight, "inversepower");
		setWeightFunction([power](const double& x) { return pow(x, - power); });
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

NestedSampler::NestedSampler() {
	setType(SamplerType::Nested);
	logEvidence = - std::numeric_limits<double>::infinity();
	logWeight = 0;
}

NestedSampler::NestedSampler(unsigned int nLivePoints) {
	setType(SamplerType::Nested);
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

	// find the live point with the lowest likelihood
	auto minIt = std::min_element(liveLikelihoods.begin(), liveLikelihoods.end());
	size_t minIdx = std::distance(liveLikelihoods.begin(), minIt);
	double minLikelihood = *minIt;

	// update log evidence
	double logVolume = - static_cast<double>(minIdx + 1) / nLivePoints;
	logWeight = logVolume + minLikelihood;
	logEvidence = logSumExp(logEvidence, logWeight);

	// Replace the worst live point with a new sample
	double newPoint;
	double newLikelihood;
	do {
		newPoint = random.randUniform(range.first, range.second);
		newLikelihood = likelihood(newPoint);
	} while (newLikelihood <= minLikelihood);

	livePoints[minIdx] = newPoint;
	liveLikelihoods[minIdx] = newLikelihood;

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


///////////////////////////////////////////////////////////////////////////////////////////////////

MCMCSampler::MCMCSampler() {
	setType(SamplerType::MCMC);
	update(0., 1.);
}

MCMCSampler::MCMCSampler(unsigned int nSteps, double stepSize) {
	setType(SamplerType::MCMC);
	setNumberOfSteps(nSteps);
	setStepSize(stepSize);
	update(0., 1.);
}

MCMCSampler::MCMCSampler(ref_ptr<Histogram1D> h, unsigned int nSteps, double stepSize) {
	setType(SamplerType::MCMC);
	setDistribution(h);
	setNumberOfSteps(nSteps);
	setStepSize(stepSize);
	update(0., 1.);
	createPDF();
}

void MCMCSampler::setNumberOfSteps(unsigned int n) {
	nSteps = n;
}

void MCMCSampler::setStepSize(double s) {
	stepSize = s;
}

void MCMCSampler::createPDF() {
	pdf = [this](const double& x) { return histogram->interpolateAt(x); };
}

void MCMCSampler::update(double x, double w) const {
	currentSample = x;
	currentWeight = w;
}

unsigned int MCMCSampler::getNumberOfSteps() const {
	return nSteps;
}

double MCMCSampler::getStepSize() const {
	return stepSize;
}

std::pair<double, double> MCMCSampler::getSample(Random& random, const std::pair<double, double>& range) const {
	// initialise the chain
	if (currentSample == 0) {
		double r;
		double w;
		do {
			r = random.randUniform(range.first, range.second);
			w = histogram->interpolateAt(r);
		} while (w == 0 or isnan(w));

		update(r, w);
	}

	for (unsigned int i = 0; i < nSteps; ++i) {
		// propose a new sample (random walk)
		double proposedSample = currentSample + random.randUniform(- stepSize, stepSize);

		// ensure the proposed sample is within the range
		if (proposedSample < range.first or proposedSample > range.second) {
			continue;
		}

		// compute the weight of the proposed sample
		double proposedWeight = histogram->interpolateAt(proposedSample);
		if (proposedWeight == 0 or isnan(proposedWeight) or currentWeight == 0) {
			continue;
		}

		// acceptance probability
		double acceptanceProbability = proposedWeight / currentWeight;

		// accept or reject the proposed sample
		if (random.rand() < acceptanceProbability) {
			update(proposedSample, proposedWeight);
			return std::make_pair(proposedSample, proposedWeight);
		}
	}

	cout << "Sample: " << currentSample << " Weight: " << currentWeight << endl;

	return std::make_pair(currentSample, currentWeight);
}

void MCMCSampler::reset() {
	currentSample = 0;
	currentWeight = 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

AdaptiveMCMCSampler::AdaptiveMCMCSampler() {
	setType(SamplerType::AdaptiveMCMC);
	currentSample = 0;
	currentWeight = 0;
	acceptedSamples = 0;
	acceptanceRate = 0;
}

AdaptiveMCMCSampler::AdaptiveMCMCSampler(unsigned int nSteps, double stepSize, double adaptationRate) {
	setType(SamplerType::AdaptiveMCMC);
	setNumberOfSteps(nSteps);
	setAdaptationRate(adaptationRate);
	currentSample = 0;
	currentWeight = 0;
	acceptedSamples = 0;
	acceptanceRate = 0;
}

AdaptiveMCMCSampler::AdaptiveMCMCSampler(ref_ptr<Histogram1D> h, unsigned int nSteps, double stepSize, double adaptationRate) {
	setDistribution(h);
	setType(SamplerType::AdaptiveMCMC);
	setNumberOfSteps(nSteps);
	setAdaptationRate(adaptationRate);
	currentSample = 0;
	currentWeight = 0;
	acceptedSamples = 0;
	acceptanceRate = 0;
}

void AdaptiveMCMCSampler::setNumberOfSteps(unsigned int n) {
	nSteps = n;
}

void AdaptiveMCMCSampler::setStepSize(double s) {
	stepSize = s;
}

void AdaptiveMCMCSampler::setAdaptationRate(double a) {
	adaptationRate = a;
}

void AdaptiveMCMCSampler::createPDF() {
	pdf = [this](const double& x) { return histogram->interpolateAt(x); };
}

unsigned int AdaptiveMCMCSampler::getNumberOfSteps() const {
	return nSteps;
}

double AdaptiveMCMCSampler::getStepSize() const {
	return stepSize;
}

double AdaptiveMCMCSampler::getAdaptationRate() const {
	return adaptationRate;
}

std::pair<double, double> AdaptiveMCMCSampler::getSample(Random& random, const std::pair<double, double>& range) const {
	// initialise the chain
	if (currentSample == 0) {
		currentSample = random.randUniform(range.first, range.second);
		currentWeight = pdf(currentSample);
	}

	for (unsigned int i = 0; i < nSteps; ++i) {
		// propose a new sample
		double proposedSample = currentSample + random.randUniform(- stepSize, stepSize);

		// ensure the proposed sample is within the range
		if (proposedSample < range.first or  proposedSample > range.second) {
			continue;
		}

		// compute the weight of the proposed sample
		double proposedWeight = pdf(proposedSample);

		// Acceptance probability
		double acceptanceProbability = proposedWeight / currentWeight;

		// accept or reject the proposed sample
		if (random.randUniform(0, 1) < acceptanceProbability) {
			currentSample = proposedSample;
			currentWeight = proposedWeight;
			acceptedSamples++;
		}

		// adapt the step size
		acceptanceRate = static_cast<double>(acceptedSamples) / (i + 1);
		// stepSize *= exp(adaptationRate * (acceptanceRate -  acceptanceRate0));
		stepSize *= exp(adaptationRate * (acceptanceRate - 0.234));
	}

	return std::make_pair(currentSample, currentWeight);
}

void AdaptiveMCMCSampler::reset() {
	currentSample = 0;
	currentWeight = 0;
	acceptedSamples = 0;
	acceptanceRate = 0;
}



} // namespace livpropa
