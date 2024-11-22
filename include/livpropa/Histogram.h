#ifndef LIVPROPA_HISTOGRAM_H
#define LIVPROPA_HISTOGRAM_H

#include <algorithm>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

#include "livpropa/Common.h"

namespace livpropa {



class AbstractBin1D: public Referenced {
	protected:
		double left;
		double right;
		double centre;

	public:
		virtual ~AbstractBin1D() = default;
		void setCentre(double c);
		void setEdges(double l, double r);
		std::pair<double, double> getEdges() const;
		double getLeftEdge() const;
		double getRightEdge() const;
		double getCentre() const;
		double getWidth() const;
		bool isInBin(const double& v) const;
};


class Bin1DLin : public AbstractBin1D {
	public:
		Bin1DLin(double l, double r) {
			setEdges(l, r);
			setCentre((l + r) / 2.);
		}
};


enum class LogBase {
	e,
	two,
	ten
};

inline double getLogBase(const LogBase& b) {
	switch (b) {
		case LogBase::e:
			return M_E;
		case LogBase::two:
			return 2.;
		case LogBase::ten:
			return 10.;
	}
}


template<LogBase B>
class Bin1DLog : public AbstractBin1D {
	public:
		Bin1DLog(double l, double r) {
			setEdges(l, r);

			double b = getBase();
			setCentre(pow(b, (logBase(l, b) + logBase(r, b)) / 2.));
		}

		double getBase() const { 
			return getLogBase(B);
		}
};

typedef Bin1DLog<LogBase::ten> Bin1DLog10;
typedef Bin1DLog<LogBase::two> Bin1DLog2;
typedef Bin1DLog<LogBase::e> Bin1DLogE;



///////////////////////////////////

class AbstractHistogram1D : public Referenced {
	public:
		typedef ref_ptr<AbstractBin1D> Bin;

	protected:
		vector<Bin> bins;
		vector<double> contents;
		vector<double> weights;
		unsigned int nBins;

	public:
		virtual ~AbstractHistogram1D() = default;
		bool isInRange(const double& v) const;
		size_t getBinIndex(const double& v) const;
		unsigned int getNumberOfBins() const;
		double leftEdge() const;
		double rightEdge() const;
		void setBinContent(const size_t& idx, const double& value);
		void setBinContents(const std::vector<double>& values);
		double getBinContent(const size_t& i) const;
		double getBinCentre(const size_t& i) const;
		vector<double> getBinEdges() const;
		vector<double> getBinCentres() const;
		vector<double> getBinContents() const;
		void push(const double& v, const double& w = 1);
		void fill(const std::vector<double>& v, const std::vector<double>& w = {});
		void normalise(double norm);
		double sum() const;
		double integrate() const;
		double operator[](const size_t& i) const;
		// virtual double interpolateAt(const double& v) const = 0;
};

template<class B>
class Histogram1 : public AbstractHistogram1D {
	public:
		Histogram1();
		Histogram1(double vMin, double vMax, unsigned int n);
		~Histogram1();
		void setBins(vector<Bin> bins);
		bool isInRange(const double& v) const;
		// virtual double interpolateAt(const double& v) const = 0;
};

typedef Histogram1<Bin1DLin> Histogram1DLin;
typedef Histogram1<Bin1DLog10> Histogram1DLog10;


// class Histogram1DLin : public Histogram1<Bin1DLin> {
// 	public:
// 		Histogram1DLin(double vMin, double vMax, unsigned int n);
// 		double interpolateAt(const double& v) const;
// 		// friend std::ostream& operator<<(std::ostream& os, const Histogram1DLin& h);
// };

// class Histogram1DLog10 : public Histogram1 {
// 	public:
// 		Histogram1DLog10(double vMin, double vMax, unsigned int n);
// 		double interpolateAt(const double& v) const;
// 		friend std::ostream& operator<<(std::ostream& os, const Histogram1DLog10& h);
// };



// typedef Histogram1<BinLinear1> Histogram1Dlin;
// typedef Histogram1<BinLogarithmic1> Histogram1Dlog;

// Histogram1<Bin1Linear> : public Histogram1<Bin1> {
// 	public:
// 		Histogram1(unsigned int n);
// };





class Histogram1D : public Referenced {
	protected:
		std::vector<double> edges;
		std::vector<double> centres;
		std::vector<double> contents;
		std::vector<double> cdf;
		unsigned int nBins;
		std::string scale;
		bool logCDF;
		double minLogCDF;

	public:
		Histogram1D(std::string scale = "lin", bool logCDF = false, double minLogCDF = -14.);
		Histogram1D(double vmin, double vmax, unsigned int n, std::string scale = "lin", bool logCDF = false, double minLogCDF = -14.);
		Histogram1D(const Histogram1D& h);
		Histogram1D(Histogram1D&& h) noexcept;
		~Histogram1D();
		void initBins(double vmin, double vmax, unsigned int n);
		void prepareCDF();
		void setScale(std::string s);
		void setLogCDF(bool log);
		void setMinimumLogCDF(double v);
		void setNumberOfBins(unsigned int n);
		void setBinContents(std::vector<double> values);
		std::string getScale() const;
		bool isCDFLogarithmic() const;
		double getMinimumLogCDF() const;
		std::vector<double> getBinEdges() const;
		std::vector<double> getBinWidths() const;
		std::vector<double> getBinCentres() const;
		std::vector<double> getBinContents() const;
		std::vector<double> getCDF() const;
		size_t getBinIndex(const double& v) const;
		unsigned int getNumberOfBins() const;
		std::pair<double, double> getEdgesOfBin(const size_t& i) const;
		double getBinCentre(const size_t& i) const;
		double getBinContent(const size_t& i) const;
		void setBinContent(size_t i, double v);
		void push(double v, double w = 1);
		void fill(const std::vector<double>& v, const std::vector<double>& w = {});
		void normalise(double norm);
		double sum() const;
		double integrate() const;
		double interpolateAt(const double& v) const;
		double getSample(bool binned = true) const;
		double getSample(const double& xMin, const double& xMax, bool binned = true) const;
		double getSample(const std::pair<double, double>& range, bool binned = true) const;
		double leftEdge() const;
		double rightEdge() const;
		void clear();
		void reset();
		double operator[](const size_t& i) const;
		Histogram1D& operator=(const Histogram1D& h);
		static std::vector<double> computeBinCentres(const std::vector<double>& binEdges, const std::string& scale);
		friend std::ostream& operator<<(std::ostream& os, const Histogram1D& h);

	private:
		static std::vector<double> computeBinEdgesLinear(double vmin, double vmax, unsigned int n);
		static std::vector<double> computeBinEdgesLog10(double vmin, double vmax, unsigned int n);
		static std::vector<double> computeBinCentresLinear(const std::vector<double>& binEdges, const std::string& scale);
		static std::vector<double> computeBinCentresLog10(const std::vector<double>& binEdges, const std::string& scale);
		static double getIntegral(const std::vector<double>& y, const std::vector<double>& dx);

		// double samplingInverse(Random& random) const;
		// double samplingImportance(Random& random, std::function<double(double)> weightFunction) const;
		
		// double samplingImportance() const;
};

} // namespace livpropa

#endif // LIVPROPA_HISTOGRAM_H