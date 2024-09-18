#ifndef LIVPROPA_HISTOGRAM_H
#define LIVPROPA_HISTOGRAM_H

#include <algorithm>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

#include "livpropa/Common.h"

namespace livpropa {

class Histogram1D : public Referenced {
	protected:
		std::vector<double> edges;
		std::vector<double> centres;
		std::vector<double> contents;
		std::vector<double> cdf;
		unsigned int nBins;
		std::string scale;

	public:
		Histogram1D(std::string scale = "lin");
		Histogram1D(double vmin, double vmax, unsigned int n, std::string scale = "lin");
		Histogram1D(const Histogram1D& h);
		Histogram1D(Histogram1D&& h) noexcept;
		~Histogram1D();
		void initBins(double vmin, double vmax, unsigned int n);
		void prepareCDF();
		void setScale(std::string s);
		void setNumberOfBins(unsigned int n);
		void setBinContents(std::vector<double> values);
		std::string getScale() const;
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
};

} // namespace livpropa

#endif // LIVPROPA_HISTOGRAM_H