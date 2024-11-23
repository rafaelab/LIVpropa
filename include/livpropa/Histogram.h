#ifndef LIVPROPA_HISTOGRAM_H
#define LIVPROPA_HISTOGRAM_H

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

#include "livpropa/Common.h"



namespace livpropa {



/**
 @class Bin1D
 @brief Abstract class for a 1D bin.
 */
class Bin1D: public Referenced {
	protected:
		double left;
		double right;
		double centre;

	public:
		virtual ~Bin1D() = default;
		void setCentre(double c);
		void setEdges(double l, double r);
		std::pair<double, double> getEdges() const;
		double getLeftEdge() const;
		double getRightEdge() const;
		double getCentre() const;
		double getWidth() const;
		bool isInBin(const double& v) const;
		double randUniform(Random& random = Random::instance()) const;
		virtual double rand(Random& random = Random::instance()) const = 0;
		virtual double directTransformation(const double& v) const = 0;
		virtual double inverseTransformation(const double& v) const = 0;
};

/**
 @class Bin1DLin
 @brief Linear 1D bin.
*/
class Bin1DLin : public Bin1D {
	public:
		Bin1DLin(double l, double r);
		double rand(Random& random = Random::instance()) const;
		double directTransformation(const double& v) const;
		double inverseTransformation(const double& v) const;
};


/**
 @class LogBase
 @brief Enumeration for the base of the logarithm.
 This class is only implemented to be used later as a template parameter.
 */
enum class LogBase {
	e,
	two,
	ten
};


/**
 @class Bin1DLog
 @brief Logarithmic 1D bin.
 The base of the logarithm is given as a template parameter.
 */
template<LogBase B>
class Bin1DLog : public Bin1D {
	public:
		Bin1DLog(double l, double r);
		double getBase() const;
		double rand(Random& random = Random::instance()) const;
		double directTransformation(const double& v) const;
		double inverseTransformation(const double& v) const;
};

typedef Bin1DLog<LogBase::ten> Bin1DLog10;
typedef Bin1DLog<LogBase::two> Bin1DLog2;
typedef Bin1DLog<LogBase::e> Bin1DLogE;




/**
 @class Histogram1D
 @brief Abstract class for a 1D histogram.
 This class is meant to be used as a base class for specific implementations of 1D histograms.
 Nevertheless, it already provides most functionalities that are common to all 1D histograms.
 */
class Histogram1D : public Referenced {
	public:
		typedef ref_ptr<Bin1D> Bin;

	protected:
		vector<Bin> bins;
		vector<double> contents;
		vector<double> weights;
		unsigned int nBins;

	public:
		virtual ~Histogram1D() = default;
		void setBinContent(const size_t& idx, const double& value);
		void setBinContents(const std::vector<double>& values);
		bool isInRange(const double& v) const;
		unsigned int getNumberOfBins() const;
		double leftEdge() const;
		double rightEdge() const;
		Bin getBin(const size_t& i) const;
		size_t getBinIndex(const double& v) const;
		double getBinContent(const size_t& i) const;
		double getBinCentre(const size_t& i) const;
		double getBinWidth(const size_t& i) const;
		vector<double> getBinEdges() const;
		vector<double> getBinCentres() const;
		vector<double> getBinContents() const;
		void push(const double& v, const double& w = 1);
		void fill(const std::vector<double>& v, const std::vector<double>& w = {});
		void normalise(double norm);
		double sum() const;
		double integrate() const;
		double operator[](const size_t& i) const;
		virtual double interpolateAt(const double& v) const = 0;
		void reset();
		void clear();
};

/**
 @class Histogram1
 @brief 1D histogram.
 This class is a template class that can be instantiated with different types of bins.
 It is meant to be used as a concrete implementation of a 1D histogram.
 */
template<class B>
class Histogram1 : public Histogram1D {
	public:
		using Histogram1D::Bin;

	public:
		Histogram1();
		Histogram1(double vMin, double vMax, unsigned int n);
		Histogram1(const Histogram1<B>& h);
		Histogram1(Histogram1<B>&& h) noexcept;
		~Histogram1();
		void setBins(vector<Bin> bins);
		double directTransformation(const double& v) const;
		double inverseTransformation(const double& v) const;
		double interpolateAt(const double& v) const;
};

typedef Histogram1<Bin1DLin> Histogram1DLin;
typedef Histogram1<Bin1DLog10> Histogram1DLog10;



} // namespace livpropa



#endif // LIVPROPA_HISTOGRAM_H