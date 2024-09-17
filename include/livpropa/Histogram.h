#ifndef LIVPROPA_HISTOGRAM_H
#define LIVPROPA_HISTOGRAM_H

#include <algorithm>
#include <iterator>
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
		using BinIterator = typename vector<double>::iterator;

	protected:
		vector<double> edges;
		vector<double> centres;
		vector<double> widths;
		vector<double> contents;
		vector<double> cdf;
		unsigned int nBins;
		string scale;

	public:
		Histogram1D(string scale = "lin") {
			setScale(scale);
		}

		Histogram1D(double vmin, double vmax, unsigned int n, string scale = "lin") {
			setScale(scale);
			setNumberOfBins(n);
			initBins(vmin, vmax, n);
		}

		Histogram1D(const Histogram1D& h) {
			reset();

			nBins = h.getNumberOfBins();

			for (size_t i = 0; i < nBins; i++) {
				edges.push_back(h.getBinEdges()[i]);
				centres.push_back(h.getBinCentres()[i]);
				widths.push_back(h.getBinWidths()[i]);
				contents.push_back(h.getBinContents()[i]);
			}
			edges.push_back(h.getBinEdges()[nBins]);

			edges.resize(nBins + 1);
			centres.resize(nBins);
			widths.resize(nBins);
			contents.resize(nBins);
		}

		Histogram1D(Histogram1D&& h) noexcept {
			*this = std::move(h);
		}

		~Histogram1D() {
			// reset();
		}

		void initBins(double vmin, double vmax, unsigned int n) {
			if (n <= 0)
				throw std::runtime_error("Number of bins must be greater than zero.");

			reset();

			// reserve memory to avoid multiple reallocations
			edges.reserve(n + 1);
			centres.reserve(n);
			widths.reserve(n);
			contents.reserve(n);

			if (scale == "log10") {
				initBinsLog10(vmin, vmax, nBins);
			} else {
				initBinsLinear(vmin, vmax, nBins);
			}


			edges.resize(nBins + 1);
			centres.resize(nBins);
			widths.resize(nBins);
			contents.resize(nBins, 0.);
		}

		void prepareCDF() {
			if (edges.size() < 2)
				throw std::runtime_error("Number of bin edges must be greater than one.");

			cdf.clear();
			
			// first entry of CDF is 0
			cdf.push_back(0);

			double cum = 0;
			for (size_t i = 1; i < nBins; i++) {
				cum += contents[i - 1];
				cdf.push_back(cum);
			}

			double tot = sum();
			for (size_t i = 0; i < nBins; i++) {
				cdf[i] /= tot;
			}

			// last entry of CDF is 0
			cdf.push_back(1);

		}	

		void setScale(string s) {
			if (s == "log")
				scale = "log10";
			else
				scale = s;

			if (s != "log" && s != "log10" && s != "lin")
				throw std::runtime_error("Unknown scale " + s + ".");;
		}

		void setNumberOfBins(unsigned int n) {
			nBins = n;
		}

		void setBinContents(vector<double> values) {
			if (values.size() != contents.size())
				throw runtime_error("Cannot set bin contents: number of values does not match number of bins.");

			contents = values;
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

		vector<double> getCDF() const {
			return cdf;
		}

		size_t getBinIndex(const double& v) const {
			if  (v < edges.front() or v > edges.back())
				return -1;

			auto bin = std::lower_bound(edges.begin(), edges.end(), v);
			
			// returns the lower edge of the bin (hence -1)
			return bin - edges.begin() - 1;
			// return std::distance(edges.begin(), bin) - 1;
		}

		unsigned int getNumberOfBins() const {
			return nBins;
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

		void push(double v, double w = 1) {
			size_t idx = getBinIndex(v);
			if (idx >= 0)
				contents[idx] += w; 			
		}

		void fill(const vector<double>& v, const vector<double>& w = {}) {
			if ((w.size() != 0) and (w.size() != v.size()))
				throw runtime_error("Cannot fill histogram: number of weights does not match number of values.");

			for (size_t i = 0; i < v.size(); i++) {
				push(v[i], w.size() == 0 ? 1. : w[i]);
			}
		}

		void normalise(double norm) {
			for (size_t i = 0; i < getNumberOfBins(); i++) {
				contents[i] /= norm;
			}
		}

		double sum() const {
			double sum = 0;
			for (size_t i = 0; i < getNumberOfBins(); i++) {
				sum += contents[i];
			}

			return sum;
		}

		double integrate() const {
			return getIntegral(contents, widths);
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

		double getSample() const {
			if (cdf.empty())
				throw runtime_error("Cannot get sample: CDF is empty. Call prepareCDF() first.");

			Random& random = Random::instance();
			double r = random.rand();

			auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
			size_t bin = std::distance(cdf.begin(), it) - 1;

			size_t binL = bin;
			size_t binU = bin + 1;

			return random.randUniform(edges[binL], edges[binU]);
		}

		double getSample(const double& xMin, const double& xMax) const {
			double vMin = (xMin < edges.front()) ? edges.front() : xMin;
			double vMax = (xMax > edges.back()) ? edges.back() : xMax;

			size_t idxL = getBinIndex(vMin);
			size_t idxU = getBinIndex(vMax) + 1;

			

			Histogram1D* h = new Histogram1D(edges[idxL], edges[idxU], idxU - idxL, scale);
			for (size_t i = idxL; i < idxU; ++i) {
				size_t idx = i - idxL;
				h->setBinContent(idx, contents[i]);
			}
			h->prepareCDF();

			double r = h->getSample();
			
			Random& random = crpropa::Random::instance();
			if (idxL == 0 or idxU == getNumberOfBins()) {
				if (r < vMin or r > edges[idxL]) 
					return random.randUniform(edges[idxL], vMin);
				if (r > vMax or r < edges[idxU])
					return random.randUniform(vMax, edges[idxU]);
			}
			delete h;

			return r;
		}

		void clear() {
			std::fill(contents.begin(), contents.end(), 0.);
		}

		void reset() {
			contents.clear();
			edges.clear();
			centres.clear();
			widths.clear();	
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

			// edges = std::move(h.edges);
			// centres = std::move(h.centres);
			// widths = std::move(h.widths);
			// contents = std::move(h.contents);
			// scale = std::move(h.scale);

			edges.resize(n + 1);
			centres.resize(n);
			widths.resize(n);
			contents.resize(n);
		
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

		friend std::ostream& operator<<(std::ostream& os, const Histogram1D& h) {
			os << "Histogram1D: " << endl;
			os << "  . scale: " << h.getScale() << endl;
			os << "  . number of bins: " << h.getNumberOfBins() << endl;
			os << "  . bin edges: " << h.getBinEdges().front() << ", " << h.getBinEdges().back() << endl;
			return os;
		}

	private:
		void initBinsLinear(double vmin, double vmax, unsigned int n) {
			for (size_t i = 0; i <= n; i++) {
				double v = vmin + i * (vmax - vmin) / n;
				edges.push_back(v);
			}
			for (size_t i = 0; i < n; i++) {
				centres.push_back((edges[i + 1] + edges[i]) / 2.);
				widths.push_back(edges[i + 1] - edges[i]);
			}
		}

		void initBinsLog10(double vmin, double vmax, unsigned int n) {
			vmin = log10(vmin);
			vmax = log10(vmax);
			for (size_t i = 0; i <= n; i++) {
				double v = vmin + i * (vmax - vmin) / n;
				edges.push_back(pow(10, v));
			}
			for (size_t i = 0; i < n; i++) {
				centres.push_back(pow(10, (log10(edges[i + 1]) + log10(edges[i])) / 2.));
				widths.push_back(edges[i + 1] - edges[i]);
			}
		}
		
		static double getIntegral(const vector<double>& y, const vector<double>& dx) {
			if (y.size() != dx.size())
				throw runtime_error("Cannot compute integral: number of y-values does not match number of bin widths.");

			double integral = 0;
			for (size_t i = 0; i < y.size(); i++) {
				integral += y[i] * dx[i];
			}

			return integral;
		}

		// static vector<double> getCDF(const vector<double>& y, const vector<double>& dx) {		
		// 	vector<double> cdf;
			
		// 	double cum = 0;
		// 	for (size_t i = 1; i < y.size(); i++) {
		// 		// cum += y[i] * dx[i];
		// 		cum += y[i - 1];
		// 		cdf.push_back(cum);
		// 	}

		// 	for (size_t i = 0; i < cdf.size(); i++) {
		// 		cdf[i] /= cdf.back();
		// 	}

		// 	return cdf;
		// }

		// static vector<double> getPDF(const vector<double>& y, const vector<double>& dx) {
		// 	double integral = getIntegral(y, dx);
			
		// 	vector<double> pdf;
		// 	for (size_t i = 0; i < y.size(); i++) {
		// 		pdf.push_back(y[i] / dx[i] / integral);
		// 	}

		// 	return pdf;
		// }

};





} // namespace livpropa

#endif // LIVPROPA_HISTOGRAM_H