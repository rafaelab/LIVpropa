#ifndef LIVPROPA_SAMPLER_H
#define LIVPROPA_SAMPLER_H

#include <functional>
#include <vector>
#include <stdexcept>

#include "livpropa/Common.h"
#include "livpropa/Histogram.h"



namespace livpropa {


enum class SamplerType {
	Inverse,
	Rejection,
	Importance
};


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
		virtual string getNameTag() const = 0;
		virtual std::pair<double, double> getSample(Random& random =  Random::instance(), const std::pair<double, double>& range = {0, 1}) const = 0;
		void setType(SamplerType t);
		SamplerType getType() const;
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
		string getNameTag() const;
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
		string getNameTag() const;
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
		string getNameTag() const;
		ref_ptr<Histogram1D> getProposalPDF() const;
		std::function<double(double)> getWeightFunction() const;
		void computeCDF();
		std::pair<double, double> getSample(Random& random = Random::instance(), const std::pair<double, double>& range = {0, 1}) const;
};




} // namespace livpropa

#endif // LIVPROPA_SAMPLER_H