#ifndef LIVPROPA_HISTOGRAM_H
#define LIVPROPA_HISTOGRAM_H


#include <string>
#include <vector>

#include <crpropa/Random.h>
#include <crpropa/Referenced.h>

#include "livpropa/Common.h"

namespace livpropa {




/**
 @class Histogram1D
 @brief Builds a one-dimensional histogram
 @todo Use templates for the bin and content type
 */
class Histogram1D : public Referenced {
	protected:
		std::vector<double> edges;
		std::vector<double> centres;
		std::vector<double> widths;
		std::vector<double> contents;
		std::string scale;

	public:
		Histogram1D(std::string scale = "lin") {
			setScale(scale);
		}

		Histogram1D(double vmin, double vmax, int nBins, std::string scale = "lin") {
			setScale(scale);
			initBins(vmin, vmax, nBins);
		}

		~Histogram1D() {
		}

		void initBins(double vmin, double vmax, int nBins) {
			// create empty histogram
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

		void setScale(std::string s) {
			if (s == "log")
				scale = "log10";

			if (s != "log" && s != "log10" && s != "lin")
				throw std::runtime_error("Unknown scale " + s + ".");

			scale = s;
		}

		std::string getScale() const {
			return scale;
		}

		std::vector<double> getBinEdges() const {
			return edges;
		}

		std::vector<double> getBinWidths() const {
			return widths;
		}

		std::vector<double> getBinCentres() const {
			return centres;
		}

		std::vector<double> getBinContents() const {
			return contents;
		}

		int getNumberOfBins() const {
			return contents.size();
		}

		double getSample() const {
			Random &random = Random::instance();
			size_t bin = random.randBin(contents);
			if (scale == "log10") {
				return pow(10, log10(edges[bin]) + random.rand() * log10(widths[bin]));
			} else {
				return edges[bin] + random.rand() * widths[bin];
			}
		}

		void push(double v, double w = 1) {
			std::vector<double>::const_iterator it = std::lower_bound(edges.begin(), edges.end(), v);
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

		std::pair<double, double> getBin(const size_t& i) const {
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
};


} // namespace livpropa

#endif // LIVPROPA_HISTOGRAM_H