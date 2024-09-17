#ifndef LIVPROPA_SAMPLER_H
#define LIVPROPA_SAMPLER_H

#include <algorithm>
#include <vector>

#include "livpropa/Common.h"
#include "livpropa/Histogram.h"



namespace livpropa {

 /**
 @class SamplerEvents
 @brief Abstract base class to handle individual events.
 */
class SamplerEvents: public Referenced {
	public:
		virtual double computeWeight(int id, double energy = 0, double energyFraction = 0, unsigned int counter = 0) const = 0;
};


/**
 @class SamplerEventsUniform
 @brief Throw away a fraction of the particles of a given type.
 The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 Note that the energy fraction is implemented through setters and getters, but it will constantly be modified.
 This is because the structure of computeWeight should remain constant (or variadic arguments used). 
 This will be done in the future.
 */
class SamplerEventsUniform: public SamplerEvents {
	protected:
		int particleId;
		double sampling;

	public:
		SamplerEventsUniform(int particleId, double sampling);
		void setSampling(double sampling);
		void setParticleId(int particleId);
		double getSampling() const;
		int getParticleId() const;
		double computeWeight(int id, double energy = 0, double energyFraction = 0, unsigned int counter = 0) const;
};


/**
 @class SamplerEventsEnergy
 @brief Throw away a fraction of the particles of a given type according to a given distribution.
 This distribution is defined from the derived classes.
 The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 This is essentially the implementation of the thinning of the EM interactions for index=1.
 In this case, the thinning is defined as: thinning = 1 - sampling.
 */
class SamplerEventsEnergy : public SamplerEvents {
	protected:
		int particleId;
		double sampling;

	public:
		SamplerEventsEnergy(int particleId, double sampling);
		void setSampling(double sampling);
		void setParticleId(int particleId);
		double getSampling() const;
		int getParticleId() const;
		double computeWeight(int id, double energy = 0, double energyFraction = 0, unsigned int counter = 0) const;
		virtual double weightFunction(int id, double energy = 0, double energyFraction = 0, unsigned int counter = 0) const = 0;
};


/**
 @class SamplerEventsEnergyFractionPowerLaw
 @brief Throw away a fraction of the particles of a given type according to a power law associated with the energy fraction it takes.
 The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 This is essentially the implementation of the thinning of the EM interactions for index=1.
 In this case, the thinning is defined as: thinning = 1 - sampling.

 */
class SamplerEventsEnergyFractionPowerLaw: public SamplerEventsEnergy {
	protected:
		int particleId;
		double sampling;
		double index;

	public:
		SamplerEventsEnergyFractionPowerLaw(double index, int particleId, double sampling);
		void setIndex(double index);
		double getIndex() const;
		double weightFunction(int id, double energy = 0, double energyFraction = 0, unsigned int counter = 0) const;
};


/**
 @class SamplerEventsEnergyFraction
 @brief Throw away a fraction of the particles of a given type.
 The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 This is essentially the implementation of the thinning of the EM interactions. 
 In this case, the thinning is defined as: thinning = 1 - sampling
 */
class SamplerEventsEnergyFraction: public SamplerEventsEnergyFractionPowerLaw {
	public:
		SamplerEventsEnergyFraction(int particleId, double sampling) : SamplerEventsEnergyFractionPowerLaw(1, particleId, sampling) {}
};


/**
 @class SamplerEventsNull
 @brief Useful dummy class to be used as initializer.
 */
class SamplerEventsNull : public SamplerEvents {
	public:
		SamplerEventsNull();
		double computeWeight(int id, double energy = 0, double energyFraction = 0, unsigned int counter = 0) const;
};


/**
 @class SamplerEventsList
 @brief List of objects of type "SamplerEvents".
 */
class SamplerEventsList : public SamplerEvents {
	protected:
		vector<ref_ptr<SamplerEvents>> samplers;

	public:
		SamplerEventsList();
		SamplerEventsList(vector<ref_ptr<SamplerEvents>> samplers);
		void add(SamplerEvents* samplers);
		inline void add(ref_ptr<SamplerEvents> SamplerEvents) {
			add(SamplerEvents.get());
		}
		double computeWeight(int id, double energy = 0, double energyFraction = 0, unsigned int counter = 0) const;
};


/**
 @class SamplerDistribution
 @brief Builds a histogram and only after it is built samples from it.
 */
class SamplerDistribution : public Referenced {
	public:
		virtual ~SamplerDistribution() = default;
		virtual vector<double> getSample(unsigned int nSamples) const {
			return vector<double>();
		}
		virtual unsigned int getSize() const {
			return 0;
		}
		virtual double getSample() const {
			return 0;
		}
		virtual ref_ptr<Histogram1D> getDistribution() const {
			return nullptr;
		}
		virtual double interpolateAt(const double &v) const {
			return 0;
		}
		virtual void prepareCDF() {
		}
		virtual void append(const vector<double>& v) {
		}
		virtual void push(const double& v) {
		}
		virtual void clear() {
		}
		vector<double> getSamples(unsigned int nSamples) const {
			vector<double> v;
			for (size_t i = 0; i < nSamples; i++) {
				double u = getSample();
				v.push_back(u);
			}
			return v;
		}
};


/**
 @class SamplerDistributionUniform
 @brief Samples uniformly from a given distribution.
 */
class SamplerDistributionUniform : public SamplerDistribution {
	protected:
		ref_ptr<Histogram1D> distribution;
		int datasetSize;

	public:
		SamplerDistributionUniform(double vmin, double vmax, unsigned int nBins, string scale = "lin");
		void initDistribution(double vmin, double vmax, unsigned int nBins, string scale = "lin");
		void setSize(int size);
		unsigned int getSize() const;
		void setDistribution(ref_ptr<Histogram1D> dist);
		double getSample() const;
		ref_ptr<Histogram1D> getDistribution() const;
		// vector<double> getSamples(unsigned int nSamples) const;
		double interpolateAt(const double &v) const;
		void prepareCDF();
		void append(const vector<double>& v);
		void push(const double& v);
		void clear();
};



} // namespace livpropa

#endif