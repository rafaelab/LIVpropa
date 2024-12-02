#ifndef LIVPROPA_WEIGHTER_H
#define LIVPROPA_WEIGHTER_H

#include <functional>
#include <vector>
#include <stdexcept>

#include "livpropa/Common.h"
#include "livpropa/Histogram.h"



namespace livpropa {


enum class WeighterType {
	UniformFraction,
	UniformNumber,
	List,
	Null
};


 /**
 @class Weighter
 @brief Abstract base class to handle individual events.
 */
class Weighter: public Referenced {
	protected:
		WeighterType type;

	public:
		virtual double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0, Random& random = Random::instance()) const = 0;
		virtual string getNameTag() const = 0;
		void setWeighterType(WeighterType type);
		WeighterType getType() const;
};

/**
 @class WeighterUniformFraction
 @brief Throw away a fraction of the particles of a given type.
 The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 It assumes that the total sampled energy of the secondaries (E_sampled) is a fraction of the total energy loss of the primaries.
 Therefore, the expected energy loss has to be known in advance.
 Note that the energy fraction is implemented through setters and getters, but it will constantly be modified.
 This is because the structure of computeWeight should remain constant (or variadic arguments used). 
 This will be done in the future.
 */
class WeighterUniformFraction: public Weighter {
	protected:
		int particleId;
		double samplingFraction;

	public:
		WeighterUniformFraction(int particleId, double sampling);
		void setSamplingFraction(double sampling);
		void setParticleId(int particleId);
		double getSamplingFraction() const;
		int getParticleId() const;
		double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0, Random& random = Random::instance()) const;
		string getNameTag() const;
};


/**
 @class WeighterUniformNumber
 @brief Throw away a fraction of the particles of a given type.
The sampler counts the number of events up to the desired number, and throws away the rest.
 */
class WeighterUniformNumber: public Weighter {
	protected:
		int particleId;
		unsigned int nEvents;

	public:
		WeighterUniformNumber(int particleId, unsigned int nEvents);
		void setNumberOfEvents(unsigned int nEvents);
		void setParticleId(int particleId);
		unsigned int getNumberOfEvents() const;
		int getParticleId() const;
		double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0, Random& random = Random::instance()) const;
		string getNameTag() const;
};


/**
 @class WeighterList
 @brief List of objects of type `Weighter`.
 */
class WeighterList : public Weighter {
	protected:
		std::vector<ref_ptr<Weighter>> weighters;
	
	public:
		WeighterList();
		WeighterList(std::vector<ref_ptr<Weighter>> weighters);
		void add(Weighter* weighter);
		// inline void add(ref_ptr<Weights> Weights) {
		// 	add(Weights.get());
		// }
		double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0, Random& random = Random::instance()) const;
		string getNameTag() const;
};



/**
 @class WeighterNull
 @brief Dummy object to ensure compatibility with the rest of the code.
 Always returns a weight of 1.
 */
class WeighterNull : public Weighter {
	protected:
		int particleId;
		unsigned int nEvents;

	public:
		WeighterNull();
		double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0, Random& random = Random::instance()) const;
		string getNameTag() const;
};




} // namespace livpropa

#endif // LIVPROPA_WEIGHTER_H