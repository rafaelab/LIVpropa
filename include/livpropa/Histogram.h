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
		virtual bool isLinear() const = 0;
		virtual bool isLog10() const = 0;
		virtual bool isLog2() const = 0;
		virtual bool isLn() const = 0;
};

/**
 @class Bin1DLin
 @brief Linear 1D bin.
*/
class Bin1DLin : public Bin1D {
	public:
		Bin1DLin();
		Bin1DLin(double l, double r);
		double rand(Random& random = Random::instance()) const;
		double directTransformation(const double& v) const;
		double inverseTransformation(const double& v) const;
		bool isLinear() const;
		bool isLog10() const;
		bool isLog2() const;
		bool isLn() const;
};


/**
 @class Bin1DLog
 @brief Logarithmic 1D bin.
 The base of the logarithm is given as a template parameter.
 */
class Bin1DLogarithmic : public Bin1D {
	protected:
		double base;

	public:
		Bin1DLogarithmic(double base = 10.);
		Bin1DLogarithmic(double l, double r, double base = 10);
		void setBase(double b);
		double getBase() const;
		double rand(Random& random = Random::instance()) const;
		double directTransformation(const double& v) const;
		double inverseTransformation(const double& v) const;
		bool isLinear() const;
		bool isLog10() const;
		bool isLog2() const;
		bool isLn() const;
};


class Bin1DLog2 : public Bin1DLogarithmic {
	public:
		Bin1DLog2();
		Bin1DLog2(double l, double r);
};


class Bin1DLog10 : public Bin1DLogarithmic {
	public:
		Bin1DLog10();
		Bin1DLog10(double l, double r);
};


class Bin1DLn : public Bin1DLogarithmic {
	public:
		Bin1DLn();
		Bin1DLn(double l, double r);
};


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
		mutable bool isPDF;
		mutable bool isCDF;

	public:
		virtual ~Histogram1D() = default;
		void setBinContent(const size_t& idx, const double& value);
		void setBinContents(const std::vector<double>& values);
		void setIsPDF(bool b);
		void setIsCDF(bool b);
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
		vector<Bin> getBins() const;
		bool getIsPDF() const;
		bool getIsCDF() const;
		bool isInRange(const double& v) const;
		void push(const double& v, const double& w = 1);
		void fill(const std::vector<double>& v, const std::vector<double>& w = {});
		void normalise(double norm);
		double sum() const;
		double integrate() const;
		vector<double> computeVectorCDF() const;
		void makePDF();
		void makeCDF();
		void reset();
		void clear();
		bool isLinear() const;
		bool isLog10() const;
		bool isLog2() const;
		bool isLn() const;
		bool isIrregular() const;
		virtual bool isRegular() const = 0;
		virtual double interpolateAt(const double& v) const = 0;
		virtual std::function<double(double)> getInterpolator(const std::pair<double, double>& range = {0., 1.}) const = 0;
		double operator[](const size_t& i) const;
		// friend std::ostream& operator<<(std::ostream& os, const Histogram1D& h);
};



/**
 @class Histogram1
 @brief 1D histogram.
 This class is a template class that can be instantiated with different types of bins.
 It is meant to be used as a concrete implementation of a 1D histogram.
 */
class RegularHistogram1D : public Histogram1D {
	public:
		using Histogram1D::Bin;

	public:
		virtual ~RegularHistogram1D() = default;
		void setBins(vector<Bin> bins);
		bool isRegular() const;
		double directTransformation(const double& v) const;
		double inverseTransformation(const double& v) const;
		double interpolateAt(const double& v) const;
		std::function<double(double)> getInterpolator(const std::pair<double, double>& range = {0., 1.}) const;
		ref_ptr<Histogram1D> getHistogramPDF() const;
		ref_ptr<Histogram1D> getHistogramCDF() const;
		// ref_ptr<Histogram1D>& operator=(const RegularHistogram1D& h);
		virtual ref_ptr<Histogram1D> clone() const = 0;
		friend std::ostream& operator<<(std::ostream& os, const RegularHistogram1D& h);
};


class Histogram1DLin : public RegularHistogram1D {
	public:
		Histogram1DLin();
		Histogram1DLin(double vMin, double vMax, unsigned int n);
		ref_ptr<Histogram1D> clone() const;
};

class Histogram1DLog10 : public RegularHistogram1D {
	public:
		Histogram1DLog10();
		Histogram1DLog10(double vMin, double vMax, unsigned int n);
		ref_ptr<Histogram1D> clone() const;
};

class Histogram1DLog2 : public RegularHistogram1D {
	public:
		Histogram1DLog2();
		Histogram1DLog2(double vMin, double vMax, unsigned int n);
		ref_ptr<Histogram1D> clone() const;
};

class Histogram1DLn : public RegularHistogram1D {
	public:
		Histogram1DLn();
		Histogram1DLn(double vMin, double vMax, unsigned int n);
		ref_ptr<Histogram1D> clone() const;
};


} // namespace livpropa


#endif // LIVPROPA_HISTOGRAM_H