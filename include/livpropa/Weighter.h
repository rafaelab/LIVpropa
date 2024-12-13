#ifndef LIVPROPA_WEIGHTER_H
#define LIVPROPA_WEIGHTER_H

#include <functional>
#include <vector>
#include <stdexcept>

#include <kiss/logger.h>

#include "livpropa/Common.h"
#include "livpropa/Histogram.h"
#include "livpropa/Sampler.h"



namespace livpropa {



 /**
 @class RuntimeWeighter
 @brief Abstract base class to handle individual events.
 Weights are assigned at runtime, during the simulation, in the `process` function.
 */
class RuntimeWeighter: public Referenced {
	public: 
		enum class Type {
			Null, 
			List,
			EnergyFraction,
			EnergyFractionPowerLaw,
			EnergyFractionUniform
		};

	protected:
		Type type;

	public:
		virtual ~RuntimeWeighter() = default;
		virtual double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0, Random& random = Random::instance()) const = 0;
		void setType(Type type);
		Type getType() const;
		string getNameTag() const;
		void reset();
};

/**
 @class WeighterNull
 @brief Dummy object to ensure compatibility with the rest of the code.
 Always returns a weight of 1.
 */
class WeighterNull : public RuntimeWeighter {
	protected:
		int particleId;
		unsigned int nEvents;

	public:
		WeighterNull();
		double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0, Random& random = Random::instance()) const;
};


/**
 @class WeighterList
 @brief List of objects of type `RuntimeWeighter`.
 */
class WeighterList : public RuntimeWeighter {
	protected:
		std::vector<ref_ptr<RuntimeWeighter>> weighters;
	
	public:
		WeighterList();
		WeighterList(std::vector<ref_ptr<RuntimeWeighter>> weighters);
		void add(RuntimeWeighter* weighter);
		double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0, Random& random = Random::instance()) const;
};


/**
 @class WeighterEnergyFraction
 @brief Throw away a fraction of the particles of a given type.
 The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 It assumes that the total sampled energy of the secondaries (E_sampled) is a fraction of the total energy loss of the primaries.
 Therefore, the expected energy loss has to be known in advance.
 Note that the energy fraction is implemented through setters and getters, but it will constantly be modified.
 This is because the structure of computeWeight should remain constant (or variadic arguments used). 
 This will be done in the future.
 */
class WeighterEnergyFraction : public RuntimeWeighter {
	protected:
		int particleId;
		std::function<double(double)> weightFunction;

	public:
		WeighterEnergyFraction();
		WeighterEnergyFraction(int particleId);
		WeighterEnergyFraction(int particleId, std::function<double(double)> func);
		~WeighterEnergyFraction();
		void setParticleId(int particleId);
		void setWeightFunction(std::function<double(double)> func);
		int getParticleId() const;
		std::function<double(double)> getWeightFunction() const;
		double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0, Random& random = Random::instance()) const;

};


/**
 @class WeighterEnergyFractionUniform
 */
class WeighterEnergyFractionUniform : public WeighterEnergyFraction {
	protected:
		double samplingFraction;

	public:
		WeighterEnergyFractionUniform();
		WeighterEnergyFractionUniform(int particleId, double sampling);
		void setSamplingFraction(double sampling);
		double getSamplingFraction() const;
		using RuntimeWeighter::setType;
		using RuntimeWeighter::getType;
};

/**
 @class WeighterEnergyFractionPowerLaw
 */
class WeighterEnergyFractionPowerLaw : public WeighterEnergyFraction {
	protected:
		double exponent;
		
	public:
		WeighterEnergyFractionPowerLaw();
		WeighterEnergyFractionPowerLaw(int particleId, double s);
		void setExponent(double v);
		double getExponent() const;
		using RuntimeWeighter::setType;
		using RuntimeWeighter::getType;
};





class PosterioriWeighterRandomEvents;
class PosterioriWeighterHistogramBins;

 /**
 @class PosterioriWeighter
 @brief Abstract base class to handle "a posteriori" weights.
 These weights are useful when the weight of an event is not known in advance, such that the whole computation must be done anyways, and only a sample of the secondary events is needed at the ended.
 Weights are assigned at runtime, during the simulation, in the `process` function.
 */
class PosterioriWeighter: public Referenced {
	public: 
		enum class Type {
			RandomEvents,
			HistogramBins
		};

	protected:
		Type type;

	public:
		virtual ~PosterioriWeighter() = default;
		void setType(Type type);
		Type getType() const;
		string getNameTag() const;
		virtual void push(const double& v, const double& w) = 0;
		virtual vector<std::pair<double, double>> getEvents(Random& random = Random::instance()) const = 0;
		virtual void reset() = 0;
};


/**
 @class PosterioriWeighterRandomEvents
 @brief Randomly select the events to be stored, according to a given sampling algorithm.
 Note that, if some kind of sampling was used to select the events to be emitted, there might be some convergence issues.
 */
class PosterioriWeighterRandomEvents: public PosterioriWeighter {
	protected:
		unsigned int nEvents;
		mutable ref_ptr<Sampler> sampler;
		mutable unsigned int nEntries;
		mutable double sumWeights;

	public:
		PosterioriWeighterRandomEvents(unsigned int nEvents, ref_ptr<Sampler> sampler = nullptr);
		void setNumberOfEvents(unsigned int nEvents);
		void setSampler(ref_ptr<Sampler> s);
		void push(const double& v, const double& w);
		unsigned int getNumberOfEvents() const;
		ref_ptr<Sampler> getSampler() const;
		ref_ptr<Histogram1D> getHistogram() const;
		double computeNormalisation() const;
		std::pair<double, double> getEvent(Random& random = Random::instance()) const;
		vector<std::pair<double, double>> getEvents(Random& random = Random::instance()) const;
		void reset();
};



/**
 @class PosterioriWeighterHistogramBins
 @brief Stores the full distribution of events.
 The events are then chose deterministically, taking the centre of each bin.
 */
class PosterioriWeighterHistogramBins: public PosterioriWeighter {
	protected:
		mutable ref_ptr<Histogram1D> histogram;
		mutable unsigned int nEntries;
		mutable double sumWeights;

	public:
		PosterioriWeighterHistogramBins(ref_ptr<Histogram1D> histogram = nullptr);
		void setHistogram(ref_ptr<Histogram1D> h);
		void push(const double& v, const double& w);
		ref_ptr<Histogram1D> getHistogram() const;
		double computeNormalisation() const;
		vector<std::pair<double, double>> getEvents(Random& random = Random::instance()) const;
		void reset();
};




} // namespace livpropa

#endif // LIVPROPA_WEIGHTER_H