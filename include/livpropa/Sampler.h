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
		ref_ptr<AbstractHistogram1D> histogram;
		vector<double> cdf;

	public:
		virtual ~Sampler() = default;
		virtual double getSample(Random& random) const = 0;
		// virtual vector<double> getSamples(unsigned int nSamples) const = 0;
		void setHistogram(ref_ptr<AbstractHistogram1D> histogram);
		ref_ptr<AbstractHistogram1D> getHistogram() const;
		void computeCDF();

		// void setBinEdges(const vector<double>& edges);
		// void setCDF(const vector<double>& cdf);
		// void setScale(const string& scale);
};


/**
 @class InverseSampler
 @brief Sample from a distribution using the inverse transform method.
*/
class SamplerInverse: public Sampler {
	public:
		SamplerInverse();
		SamplerInverse(ref_ptr<AbstractHistogram1D> h);
		~SamplerInverse();
		double getSample(Random& random) const;
		// vector<double> getSamples(unsigned int nSamples) const override;
};


// // /**
// //  @class ImportanceSampler
// //  @brief Sample from a distribution using the importance sampling method.
// // */
// // class ImportanceSampler: public Sampler {
// // 	private:
// // 		vector<double> edges;
// // 		vector<double> cdf;
// // 		vector<double> proposal;
// // 		vector<double> contents;
// // 		string scale;
// // 		double total;

// // 	public:
// // 		ImportanceSampler(const vector<double>& edges, const vector<double>& cdf, const vector<double>& proposalCDF, const vector<double>& contents, double total, const string& scale);
// // 		double sample() const override;
// // };




} // namespace livpropa

#endif // LIVPROPA_SAMPLER_H