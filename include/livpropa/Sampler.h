#ifndef LIVPROPA_SAMPLER_H
#define LIVPROPA_SAMPLER_H

#include <functional>
#include <unordered_map>
#include <vector>
#include <stdexcept>

#include "livpropa/Common.h"
#include "livpropa/Histogram.h"



namespace livpropa {


enum class SamplerType {
	Inverse,
	Rejection,
	Importance,
	Nested,
	MCMC,
	AdaptiveMCMC
};

string getSamplerNameTag(SamplerType t);


/**
 @class Sampler
 @brief Interface for sampling from a distribution.
*/
class Sampler : public Referenced {
	protected:
		ref_ptr<Histogram1D> histogram;
		vector<double> cdf;
		SamplerType type;

	public:
		virtual ~Sampler() = default;
		virtual std::pair<double, double> getSample(Random& random =  Random::instance(), const std::pair<double, double>& range = {0, 1}) const = 0;
		void setType(SamplerType t);
		SamplerType getType() const;
		string getNameTag() const;
		vector<std::pair<double, double>> getSamples(unsigned int nSamples, Random& random = Random::instance(), const std::pair<double, double>& range = {0, 1}) const;
		void setDistribution(ref_ptr<Histogram1D> histogram);
		ref_ptr<Histogram1D> getDistribution() const;
		ref_ptr<Histogram1D> getCumulativeDistribution() const;
		void computeCDF();

};


/**
 @class InverseSampler
 @brief Sample from a distribution using the inverse transform method.
*/
class InverseSampler: public Sampler {
	public:
		InverseSampler();
		InverseSampler(ref_ptr<Histogram1D> h);
		std::pair<double, double> getSample(Random& random = Random::instance(), const std::pair<double, double>& range = {0, 1}) const;
};



/**
 @class RejectionSampler
 @brief Sample from a distribution using the rejection sampling method.
 (UNTESTED)
*/
class RejectionSampler: public Sampler {
	protected:
		ref_ptr<Histogram1D> proposalPDF;
		InverseSampler inverseSampler;
		double maxRatio;

	public:
		RejectionSampler();
		RejectionSampler(ref_ptr<Histogram1D> pdf, ref_ptr<Histogram1D> proposalPDF, double maxRatio);
		void setProposalPDF(ref_ptr<Histogram1D> proposalPDF);
		ref_ptr<Histogram1D> getProposalPDF() const;
		void computeCDF();
		std::pair<double, double> getSample(Random& random = Random::instance(), const std::pair<double, double>& range = {0, 1}) const;
};


/**
 @class ImportanceSampler
 @brief Sample from a distribution using the importance sampling method.
 The default CDF function is changed here to incorporate the proposal distribution.
*/
class ImportanceSampler: public Sampler {
	protected:
		std::function<double(double)> weightFunction;
		ref_ptr<Histogram1D> proposalPDF;
		InverseSampler inverseSampler;

	public:
		ImportanceSampler();
		ImportanceSampler(std::function<double(double)> weight);
		ImportanceSampler(string weight);
		ImportanceSampler(ref_ptr<Histogram1D> pdf, ref_ptr<Histogram1D> proposalPDF, std::function<double(double)> weight);
		ImportanceSampler(ref_ptr<Sampler> sampler);
		void setProposalPDF(ref_ptr<Histogram1D> proposalPDF);
		void setWeightFunction(std::function<double(double)> func);
		ref_ptr<Histogram1D> getProposalPDF() const;
		std::function<double(double)> getWeightFunction() const;
		void computeCDF();
		std::pair<double, double> getSample(Random& random = Random::instance(), const std::pair<double, double>& range = {0, 1}) const;
	
	private:
		static double parseWeightFunctionName(const string& str, const string& pattern);
};


/**
 @class NestedSampler
 @brief Sample from a distribution using the Nested Sampling method.
*/
class NestedSampler : public Sampler {
	protected:
		unsigned int nLivePoints;
		std::function<double(double)> likelihood;
		mutable std::vector<double> livePoints;
		mutable std::vector<double> liveLikelihoods;
		mutable double logEvidence;
		mutable double logWeight;

	public:
		NestedSampler();
		NestedSampler(unsigned int nLivePoints);
		NestedSampler(ref_ptr<Histogram1D> h, unsigned int nLivePoints);
		void setNumberOfLivePoints(unsigned int n);
		void setLikelihoodFunction(std::function<double(double)> func);
		void setDistribution(ref_ptr<Histogram1D> h);
		unsigned int getNumberOfLivePoints() const;
		std::function<double(double)> getLikelihoodFunction() const;
		std::pair<double, double> getSample(Random& random = Random::instance(), const std::pair<double, double>& range = {0, 1}) const override;
		double getLogEvidence() const;

	private:
		double logSumExp(double a, double b) const;
};


/**
 @class MCMCSampler
 @brief Sample from a distribution using the Markov Chain Monte Carlo (MCMC) method.
 [UNTESTED]
*/
class MCMCSampler : public Sampler {
	protected:
		unsigned int nSteps;
		double stepSize;
		std::function<double(double)> pdf;
			InverseSampler inverseSampler;
		mutable double currentSample;
		mutable double currentWeight;

	public:
		MCMCSampler();
		MCMCSampler(unsigned int nSteps, double stepSize);
		MCMCSampler(ref_ptr<Histogram1D> h, unsigned int nSteps, double stepSize);
		void setNumberOfSteps(unsigned int n);
		void setStepSize(double s);
		void createPDF();
		void update(double x, double w) const;
		unsigned int getNumberOfSteps() const;
		double getStepSize() const;
		std::pair<double, double> getSample(Random& random = Random::instance(), const std::pair<double, double>& range = {0, 1}) const;
		void reset();
};


/**
 @class AdaptativeMCMCSampler
 @brief Sample from a distribution using the Markov Chain Monte Carlo (MCMC) method.
 [UNTESTED]
*/
class AdaptiveMCMCSampler : public Sampler {
	protected:
		unsigned int nSteps;
		double adaptationRate;
		std::function<double(double)> pdf;
		mutable double stepSize;
		mutable double currentSample;
		mutable double currentWeight;
		mutable unsigned int acceptedSamples;
		mutable double acceptanceRate;
		// constexpr static double acceptanceRate0 = 0.234;

	public:
		AdaptiveMCMCSampler();
		AdaptiveMCMCSampler(unsigned int nSteps, double stepSize, double adaptationRate);
		AdaptiveMCMCSampler(ref_ptr<Histogram1D> h, unsigned int nSteps, double stepSize, double adaptationRate);
		void setNumberOfSteps(unsigned int n);
		void setStepSize(double s);
		void setAdaptationRate(double a);
		void createPDF();
		unsigned int getNumberOfSteps() const;
		double getStepSize() const;
		double getAdaptationRate() const;
		std::pair<double, double> getSample(Random& random = Random::instance(), const std::pair<double, double>& range = {0, 1}) const;
};

} // namespace livpropa

#endif // LIVPROPA_SAMPLER_H