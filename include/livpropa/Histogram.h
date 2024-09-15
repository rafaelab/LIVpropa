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

		Histogram1D(double vmin, double vmax, int nBins, string scale = "lin") {
			setScale(scale);
			initBins(vmin, vmax, nBins);
		}

		Histogram1D(const Histogram1D& h) {
			if (getNumberOfBins() != h.getNumberOfBins())
				throw std::runtime_error("Histograms have different number of bins.");
			
			edges.clear();
			centres.clear();
			widths.clear();
			contents.clear();
			
			for (size_t i = 0; i < h.getNumberOfBins(); i++) {
				edges.push_back(h.getBinEdges()[i]);
				centres.push_back(h.getBinCentres()[i]);
				widths.push_back(h.getBinWidths()[i]);
				contents.push_back(h.getBinContents()[i]);
			}
			edges.push_back(h.getBinEdges().back());

		}

		~Histogram1D() {
		}

		void initBins(double vmin, double vmax, int nBins) {
			if (scale == "log10") {
				vmin = log10(vmin);
				vmax = log10(vmax);
				for (size_t i = 0; i < nBins + 1; i++) {
					double v = vmin + i * (vmax - vmin) / (nBins);
					edges.push_back(pow(10, v));
				}
				for (size_t i = 0; i < nBins; i++) {
					centres.push_back(pow(10, (log10(edges[i + 1]) + log10(edges[i])) / 2.));
				}
			} else {
				for (size_t i = 0; i < nBins + 1; i++) {
					double v = vmin + i * (vmax - vmin) / (nBins);
					edges.push_back(v);
				}
				for (size_t i = 0; i < nBins; i++) {
					centres.push_back((edges[i + 1] + edges[i]) / 2.);
				}
			}

			// bins start with no content
			for (size_t i = 0; i < nBins; i++) {
				widths.push_back(edges[i + 1] - edges[i]);
				contents.push_back(0);
			}
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

		int getNumberOfBins() const {
			return contents.size();
		}

		double getSample() const {
			Random &random = Random::instance();
			size_t bin = random.randBin(contents); // should be a CDF
			if (scale == "log10") {
				return pow(10, log10(edges[bin]) + random.rand() * log10(widths[bin]));
			} else {
				return edges[bin] + random.rand() * widths[bin];
			}
		}

		double getSampleInRange(const double& xMin, const double& xMax) const {
			BinIterator binL = whichBin(xMin, true);
			BinIterator binU = whichBin(xMax, false);

			vector<double> edgesInRange; 
			vector<double> contentsInRange;
			vector<double> widthsInRange;
			for (auto i = binL; i != binU; ++i) {
				edgesInRange.push_back(*i);
				contentsInRange.push_back(contents[i - edges.begin()]);
				widthsInRange.push_back(widths[i - edges.begin()]);
			}

			Random &random = Random::instance();
			size_t bin = random.randBin(contentsInRange);
			if (scale == "log10") {
				return pow(10, log10(edgesInRange[bin]) + random.rand() * log10(widthsInRange[bin]));
			} else {
				return edgesInRange[bin] + random.rand() * widthsInRange[bin];
			}

		}

		double getSampleInRange(const std::pair<double, double>& range) const {
			return getSampleInRange(range.first, range.second);
		}

		void push(double v, double w = 1) {
			BinIterator it = whichBin(v);
			if (it == edges.begin() || it == edges.end())
				return;

			size_t idx = it - edges.begin();
			contents[idx] += w; 
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
					integral +=  (contents[i] / (edges[i + 1] - edges[i]) * log(10) * centres[i]);
				}
			} else {
				for (size_t i = 0; i < getNumberOfBins(); i++) {
					// integral += contents[i] / (edges[i] - edges[i - 1]);
					integral += contents[i];
				}
			}

			return integral;
		}

		void transformToPDF() {
			double integral = integrate();
			for (size_t i = 1; i < getNumberOfBins(); i++) {
				// contents[i] /= (edges[i] - edges[i - 1]);
				contents[i] /= integral;
			}
		}

		void transformToCDF() {
			for (size_t i = 1; i < getNumberOfBins(); i++) {
				contents[i] += contents[i - 1];
			}
		}

		void clear() {
			for (size_t i = 0; i < getNumberOfBins(); i++) {
				contents[i] = 0;
			}
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

			edges.clear();
			centres.clear();
			widths.clear();
			contents.clear();
			
			for (size_t i = 0; i < h.getNumberOfBins(); i++) {
				edges.push_back(h.getBinEdges()[i]);
				centres.push_back(h.getBinCentres()[i]);
				widths.push_back(h.getBinWidths()[i]);
				contents.push_back(h.getBinContents()[i]);
			}
			edges.push_back(h.getBinEdges().back());
		
			return *this;
		}

	private:
		BinIterator whichBin(double v, bool lowerEdge = true) const {
			BinIterator it;
			if (lowerEdge) 
				return std::lower_bound(edges.begin(), edges.end(), v);
			else
				return std::upper_bound(edges.begin(), edges.end(), v);
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