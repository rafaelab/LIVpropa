#ifndef LIVPROPA_SAMPLER_H
#define LIVPROPA_SAMPLER_H

#include <functional>
#include <vector>
#include <stdexcept>

#include "livpropa/Common.h"
#include "livpropa/Histogram.h"



namespace livpropa {


/**
 @class Sampler
 @brief Interface for sampling from a distribution.
*/
class Sampler : public crpropa::Referenced {
	protected:
		ref_ptr<Histogram1D> histogram;
		vector<double> cdf;

	public:
		virtual ~Sampler() = default;
		virtual std::pair<double, double> getSample(Random& random =  Random::instance(), const std::pair<double, double>& range = {0, 1}) const = 0;
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
 @class ImportanceSampler
 @brief Sample from a distribution using the importance sampling method.
 The default CDF function is changed here to incorporate the proposal distribution.
*/
class ImportanceSampler: public Sampler {
	protected:
		ref_ptr<Histogram1D> proposalPDF;
		InverseSampler inverseSampler;

	public:
		ImportanceSampler();
		ImportanceSampler(ref_ptr<Histogram1D> pdf, ref_ptr<Histogram1D> proposalPDF);
		void setProposalPDF(ref_ptr<Histogram1D> proposalPDF);
		ref_ptr<Histogram1D> getProposalPDF() const;
		void computeCDF();
		std::pair<double, double> getSample(Random& random = Random::instance(), const std::pair<double, double>& range = {0, 1}) const;
};




} // namespace livpropa

#endif // LIVPROPA_SAMPLER_H