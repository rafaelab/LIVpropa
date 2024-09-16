#ifndef LIVPROPA_HISTOGRAM_H
#define LIVPROPA_HISTOGRAM_H

#include <algorithm>
#include <string>
#include <vector>

#include "livpropa/Common.h"



namespace livpropa {




/**
 @class Histogram1D
 @brief Builds a one-dimensional histogram
 @todo Use templates for the bin and content type
 */
class Histogram1D : public Referenced {
	private:
		using BinIterator = typename vector<double>::const_iterator;

	protected:
		vector<double> edges;
		vector<double> centres;
		vector<double> widths;
		vector<double> contents;
		string scale;

	public:
		Histogram1D(string scale = "lin") {
			setScale(scale);
		}

		Histogram1D(double vmin, double vmax, unsigned int nBins, string scale = "lin") {
			setScale(scale);
			initBins(vmin, vmax, nBins);
		}

		Histogram1D(const Histogram1D& h) {
			reset();

			unsigned int n = h.getNumberOfBins();

			for (size_t i = 0; i < n; i++) {
				edges.push_back(h.getBinEdges()[i]);
				centres.push_back(h.getBinCentres()[i]);
				widths.push_back(h.getBinWidths()[i]);
				contents.push_back(h.getBinContents()[i]);
			}
			edges.push_back(h.getBinEdges()[n]);

			edges.resize(n + 1);
			centres.resize(n);
			widths.resize(n);
			contents.resize(n);
		}

		~Histogram1D() {
		}

		void initBins(double vmin, double vmax, unsigned int nBins) {
			if (nBins <= 0)
				throw std::runtime_error("Number of bins must be greater than zero.");

			reset();

			// reserve memory to avoid multiple reallocations
			edges.reserve(nBins + 1);
			centres.reserve(nBins);
			widths.reserve(nBins);
			contents.reserve(nBins);

			if (scale == "log10") {
				initBinsLog10(vmin, vmax, nBins);
			} else {
				initBinsLinear(vmin, vmax, nBins);
			}

			// bins starts with no content
			for (size_t i = 0; i < nBins; i++) {
				contents.push_back(0.);
			}

			edges.resize(nBins + 1);
			centres.resize(nBins);
			widths.resize(nBins);
			contents.resize(nBins);
		}

		void setScale(string s) {
			if (s == "log")
				scale = "log10";

			if (s != "log" && s != "log10" && s != "lin")
				throw std::runtime_error("Unknown scale " + s + ".");

			scale = s;
		}

		string getScale() const {
			return scale;
		}

		vector<double> getBinEdges() const {
			return edges;
		}

		vector<double> getBinWidths() const {
			return widths;
		}

		vector<double> getBinCentres() const {
			return centres;
		}

		vector<double> getBinContents() const {
			return contents;
		}

		size_t getBinIndex(const double& v) const {
			BinIterator bin = whichBin(v);
			return std::distance(edges.begin(), bin);
		}

		unsigned int getNumberOfBins() const {
			return contents.size();
		}

		double drawRandom() const {
			return getRandomValue(centres, contents, widths, scale);
		}

		double drawRandomInRange(const double& xMin, const double& xMax) const {
			BinIterator binL = whichBin(xMin, true);
			BinIterator binU = whichBin(xMax, false);

			vector<double> edgesInRange;
			vector<double> centresInRange; 
			vector<double> contentsInRange;
			vector<double> widthsInRange;
			for (auto i = binL; i != binU; ++i) {
				edgesInRange.push_back(*i);
				centresInRange.push_back(centres[i - edges.begin()]);
				contentsInRange.push_back(contents[i - edges.begin()]);
				widthsInRange.push_back(widths[i - edges.begin()]);
			}

			return getRandomValue(centresInRange, contentsInRange, widthsInRange, scale);
		}

		double getSample() const {
			double v = drawRandom();
			return interpolate(v, centres, contents);
		}

		double getSampleInRange(const double& xMin, const double& xMax) const {
			double v = drawRandomInRange(xMin, xMax);
			return interpolate(v, centres, contents);
		}

		void push(double v, double w = 1) {
			BinIterator it = whichBin(v);
			if (it == edges.begin() || it == edges.end())
				return;

			size_t idx = it - edges.begin();
			contents[idx] += w; 
		}

		void fill(const vector<double>& v, const vector<double>& w = {}) {
			if (w.size() != 0) {
				if (w.size() != v.size()) {
					throw runtime_error("Cannot fill histogram: number of weights does not match number of values.");
				}
			} 
			
			for (size_t i = 0; i < v.size(); i++) {
				push(v[i], w.size() == 0 ? 1 : w[i]);
			}
			
		}

		void normalise(double norm) {
			for (size_t i = 0; i < getNumberOfBins(); i++) {
				contents[i] /= norm;
			}
		}

		double sum() {
			double sum = 0;
			for (size_t i = 0; i < getNumberOfBins(); i++) {
				sum += contents[i];
			}

			return sum;
		}

		double integrate() {
			double integral = 0;
			if (scale == "log") {
				for (size_t i = 0; i < getNumberOfBins(); i++) {
					integral += contents[i] * (log10(edges[i + 1]) - log10(edges[i]));
				}
			} else {
				for (size_t i = 0; i < getNumberOfBins(); i++) {
					integral += contents[i] * widths[i];
				}
			}

			return integral;
		}

		double interpolateAt(const double& v) const {
			if (scale == "log10") {
				vector<double> x;
				vector<double> y;
				for (size_t i = 0; i < getNumberOfBins(); i++) {
					x.push_back(log10(centres[i]));
					y.push_back(contents[i]);
				}
				return interpolate(log10(v), x, y);
			} else {
				return interpolate(v, centres, contents);
			}
		}

		void transformToPDF() {
			double totalSum = sum();
			if (totalSum == 0) 
				throw runtime_error("Cannot transform to PDF: total sum of bin contents is zero.");
			
			for (size_t i = 0; i < getNumberOfBins(); i++) {
				contents[i] /= totalSum;
			}
		}

		void transformToCDF() {
			transformToPDF();
			for (size_t i = 1; i < getNumberOfBins(); i++) {
				contents[i] += contents[i - 1];
			}
		}

		void clear() {
			for (size_t i = 0; i < getNumberOfBins(); i++) {
				contents[i] = 0;
			}
		}

		void reset() {
			contents.clear();
			edges.clear();
			centres.clear();
			widths.clear();	
		}

		void setBinContent(size_t i, double v) {
			contents[i] = v;
		}

		std::pair<double, double> getEdgesOfBin(const size_t& i) const {
			return std::make_pair(edges[i], edges[i + 1]);
		}

		double getBinCentre(const size_t& i) const {
			if (scale == "log10") {
				return pow(10, (log10(edges[i]) + log10(edges[i + 1])) / 2.);
			} else {
				return (edges[i] + edges[i + 1]) / 2.;
			}
		}

		double getBinContent(const size_t& i) const {
			return contents[i];
		}

		double operator[](const size_t& i) const {
			return getBinContent(i);
		}

		Histogram1D& operator=(const Histogram1D& h) {
			if (this == &h)
				return *this;

			reset();
			
			unsigned int n = h.getNumberOfBins();
			for (size_t i = 0; i < n; i++) {
				edges.push_back(h.getBinEdges()[i]);
				centres.push_back(h.getBinCentres()[i]);
				widths.push_back(h.getBinWidths()[i]);
				contents.push_back(h.getBinContents()[i]);
			}
			edges.push_back(h.getBinEdges()[n]);
		
			return *this;
		}

		static vector<double> computeBinCentres(const vector<double>& binEdges, const string& scale) {
			vector<double> binCentres;
			for (size_t i = 0; i < binEdges.size() - 1; i++) {
				if (scale == "log10") {
					binCentres.push_back(pow(10, (log10(binEdges[i]) + log10(binEdges[i + 1])) / 2.));
				} else {
					binCentres.push_back((binEdges[i] + binEdges[i + 1]) / 2.);
				}
			}
			return binCentres;
		}

	private:
		void initBinsLinear(double vmin, double vmax, unsigned int nBins) {
			for (size_t i = 0; i <= nBins; i++) {
				double v = vmin + i * (vmax - vmin) / nBins;
				edges.push_back(v);
			}
			for (size_t i = 0; i < nBins; i++) {
				centres.push_back((edges[i + 1] + edges[i]) / 2.);
				widths.push_back(edges[i + 1] - edges[i]);
			}
		}

		void initBinsLog10(double vmin, double vmax, unsigned int nBins) {
			vmin = log10(vmin);
			vmax = log10(vmax);
			for (size_t i = 0; i <= nBins; i++) {
				double v = vmin + i * (vmax - vmin) / nBins;
				edges.push_back(pow(10, v));
			}
			for (size_t i = 0; i < nBins; i++) {
				centres.push_back(pow(10, (log10(edges[i + 1]) + log10(edges[i])) / 2.));
				widths.push_back(edges[i + 1] - edges[i]);
			}
		}

		BinIterator whichBin(double v, bool lowerEdge = true) const {
			BinIterator it;
			if (lowerEdge) 
				return std::lower_bound(edges.begin(), edges.end(), v);
			else
				return std::upper_bound(edges.begin(), edges.end(), v);
		}
		
		static double getRandomValue(const vector<double>& binCentres, const vector<double>& binContents, const vector<double>& binWidths, const string& scale) {	
			Random& random = Random::instance();
			return interpolate(random.rand(), binContents, binCentres);
		}
};


// /** 
// This function prints Histogram-type objects in a nice way
// */
// std::ostream& operator<<(std::ostream& os, const Histogram1D& h) {
// 	os << "Histogram1D: " << endl;
// 	os << "  . scale: " << h.getScale() << endl;
// 	os << "  . number of bins: " << h.getNumberOfBins() << endl;
// 	os << "  . bin edges: " << h.getBinEdges().front() << ", " << h.getBinEdges().back() << endl;
// 	return os;
// }


} // namespace livpropa

#endif // LIVPROPA_HISTOGRAM_H