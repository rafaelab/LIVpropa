#include "livpropa/Sampler.h"

namespace livpropa {


void Sampler::setHistogram(ref_ptr<AbstractHistogram1D> h) {
	histogram = h;
}

ref_ptr<AbstractHistogram1D> Sampler::getHistogram() const {
	return histogram;
}

void Sampler::computeCDF() {

	

	vector<double> edges = histogram->getBinEdges();
	vector<double> contents = histogram->getBinContents();
	size_t nBins = histogram->getNumberOfBins();
	// string scale = histogram->getScale();
	// double minLogCDF = histogram->getMinimumLogCDF();

		if (edges.size() < 2)
		throw std::runtime_error("Number of bin edges must be greater than one.");

	// // ensure vector is empty
	// cdf.clear();

	// first entry of CDF is 0
	cdf.push_back(0);

		// if (scale == "lin") {
		// 	double cum = 0;
		// 	for (size_t i = 0; i < nBins - 1; i++) {
		// 		cum += contents[i + 1];
		// 		cdf.push_back(cum);
		// 	}
		// } else if (scale == "log10") {
		// 	double cum = 0;
		// 	for (size_t i = 1; i < nBins; i++) {
		// 		double dx = edges[i] - edges[i - 1];
		// 		double y = log10(contents[i] + pow(10, minLogCDF - 1));
		// 		// double dy = contents[i] - contents[i - 1];
		// 		// double x = centres[i - 1];
		// 		// double y = contents[i - 1];
		// 		// y = log10(y);
		// 		// double w = y;
		// 		// cum += (y * dx / 2. * w);
		// 		cum += y * dx;
		// 		cdf.push_back(cum);
		// 	}
		// }

		// double tot = histogram->sum();
		for (size_t i = 1; i < nBins; i++) {
		// 	cdf[i] /= tot;
		cdf.push_back(i / nBins);
		}

		// // do cdf->log(cdf)
		// cdf[0] = minLogCDF;
		// for (size_t i = 1; i < nBins; i++) {
		// 	cdf[i] = log10(cdf[i]);
		// }
}


// void Histogram1D::prepareCDF() {
// 	if (edges.size() < 2)
// 		throw std::runtime_error("Number of bin edges must be greater than one.");

// 	// ensure vector is empty
// 	cdf.clear();

// 	// first entry of CDF is 0
// 	cdf.push_back(0);

// 	if (logCDF) {
// 		if (scale == "lin") {
// 			double cum = 0;
// 			for (size_t i = 0; i < nBins - 1; i++) {
// 				cum += contents[i + 1];
// 				cdf.push_back(cum);
// 			}
// 		} else if (scale == "log10") {
// 			double cum = 0;
// 			for (size_t i = 1; i < nBins; i++) {
// 				double dx = edges[i] - edges[i - 1];
// 				double y = log10(contents[i] + pow(10, minLogCDF - 1));
// 				// double dy = contents[i] - contents[i - 1];
// 				// double x = centres[i - 1];
// 				// double y = contents[i - 1];
// 				// y = log10(y);
// 				// double w = y;
// 				// cum += (y * dx / 2. * w);
// 				cum += y * dx;
// 				cdf.push_back(cum);
// 			}
// 		}

// 		double tot = sum();
// 		for (size_t i = 0; i < nBins; i++) {
// 			cdf[i] /= tot;
// 		}

// 		// // do cdf->log(cdf)
// 		// cdf[0] = minLogCDF;
// 		// for (size_t i = 1; i < nBins; i++) {
// 		// 	cdf[i] = log10(cdf[i]);
// 		// }



// 	} else {
// 		if (scale == "lin") {
// 			double cum = 0;
// 			for (size_t i = 0; i < nBins - 1; i++) {
// 				cum += contents[i + 1];
// 				cdf.push_back(cum);
// 			}
// 		} else if (scale == "log10") {
// 			double cum = 0;
// 			for (size_t i = 1; i < nBins; i++) {
// 				double dx = edges[i] - edges[i - 1];
// 				double dy = contents[i] - contents[i - 1];
// 				cum += (contents[i - 1] * dx / 2.);
// 				cdf.push_back(cum);
// 			}
// 		}

// 		double tot = sum();
// 		for (size_t i = 0; i < nBins; i++) {
// 			cdf[i] /= tot;
// 		}
// 	}

// }	




SamplerInverse::SamplerInverse() {
}

SamplerInverse::SamplerInverse(ref_ptr<AbstractHistogram1D> h) {
	setHistogram(h);
	computeCDF();
}

SamplerInverse::~SamplerInverse() {
}

double SamplerInverse::getSample(Random& random) const {
	double r = random.rand();

	auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
	size_t idx = std::distance(cdf.begin(), it) - 1;

	// if (idx >= edges.size() - 1)
	// 	idx = edges.size() - 2;

	return interpolate(r, cdf, histogram->getBinEdges());
}


} // namespace livpropa
