#include "livpropa/Histogram.h"

namespace livpropa {



void AbstractBin1D::setEdges(double l, double r) {
	left = l;
	right = r;
}

void AbstractBin1D::setCentre(double c) {
	centre = c;
}

pair<double, double> AbstractBin1D::getEdges() const {
	return std::make_pair(left, right);
}

double AbstractBin1D::getLeftEdge() const {
	return left;
}

double AbstractBin1D::getRightEdge() const {
	return right;
}

double AbstractBin1D::getCentre() const {
	return centre;
}

double AbstractBin1D::getWidth() const {
	return right - left;
}

bool AbstractBin1D::isInBin(const double& v) const {
	return v >= left and v < right;
}



bool AbstractHistogram1D::isInRange(const double& v) const {
	if (v >= bins[0]->getLeftEdge() and v < bins[nBins - 1]->getRightEdge())
		return true;

	return false;
}

size_t AbstractHistogram1D::getBinIndex(const double& v) const {
	if (! isInRange(v))
		return -1;

	for (size_t i = 0; i < nBins; i++) {
		if (bins[i]->isInBin(v))
			return i;
	}

	return -1;
}

unsigned int AbstractHistogram1D::getNumberOfBins() const {
	return nBins;
}

double AbstractHistogram1D::leftEdge() const {
	return bins[0]->getLeftEdge();
}

double AbstractHistogram1D::rightEdge() const {
	return bins[nBins - 1]->getRightEdge();
}

double AbstractHistogram1D::getBinContent(const size_t& i) const {
	return contents[i];
}

double AbstractHistogram1D::getBinCentre(const size_t& i) const {
	return bins[i]->getCentre();
}

vector<double> AbstractHistogram1D::getBinEdges() const {
	vector<double> edges;
	for (size_t i = 0; i < nBins; i++) {
		edges.push_back(bins[i]->getLeftEdge());
	}
	edges.push_back(bins[nBins - 1]->getRightEdge());

	return edges;
}

vector<double> AbstractHistogram1D::getBinCentres() const {
	vector<double> centres;
	for (size_t i = 0; i < nBins; i++) {
		centres.push_back(bins[i]->getCentre());
	}

	return centres;
}

vector<double> AbstractHistogram1D::getBinContents() const {
	return contents;
}

void AbstractHistogram1D::push(const double& v, const double& w) {
	size_t i = getBinIndex(v);
	contents[i] += w;
}

void AbstractHistogram1D::fill(const vector<double>& v, const vector<double>& w) {
	if (w.size() == 0) {
		for (size_t i = 0; i < v.size(); i++) 
			push(v[i]);
	} else {
		if (v.size() != w.size())
			throw std::runtime_error("Number of values to fill bins must match number of weights.");

		for (size_t i = 0; i < v.size(); i++)
			push(v[i], w[i]);
	}
}

void AbstractHistogram1D::setBinContent(const size_t& idx, const double& value) {
	contents[idx] = value;
}

void AbstractHistogram1D::setBinContents(const vector<double>& values) {
	if (values.size() != nBins)
		throw std::runtime_error("Number of values must match number of bins.");

	contents = values;
}

void AbstractHistogram1D::normalise(double norm) {
	for (size_t i = 0; i < nBins; i++) {
		contents[i] /= norm;
	}
}

double AbstractHistogram1D::sum() const {
	double s = 0;
	for (size_t i = 0; i < nBins; i++) {
		s += contents[i];
	}

	return s;
}

double AbstractHistogram1D::integrate() const {
	double s = 0;
	for (size_t i = 0; i < nBins; i++) {
		double dx = bins[i]->getWidth();
		s += contents[i] * dx;
	}

	return s;
}

double AbstractHistogram1D::operator[](const size_t& i) const {
	return contents[i];
}




template<class B>
Histogram1<B>::Histogram1() {
	nBins = 0;
}

template<>
Histogram1<Bin1DLin>::Histogram1(double vMin, double vMax, unsigned int n) {
	nBins = n;

	vector<Bin> bs;
	double dx = (vMax - vMin) / n;
	for (size_t i = 0; i < n; i++) {
		double l = vMin + i * dx;
		double r = vMin + (i + 1) * dx;
		Bin b = new Bin1DLin(l, r);
		bs.push_back(b);
	}
	setBins(bs);
}

template<>
Histogram1<Bin1DLog10>::Histogram1(double vMin, double vMax, unsigned int n) {
	nBins = n;

	vector<Bin> bs;
	double dx = (log10(vMax / vMin)) / n;
	for (size_t i = 0; i < n; i++) {
		double l = log10(vMin) + i * dx;
		double r = log10(vMin) + (i + 1) * dx;
		l = pow(10, l);
		r = pow(10, r);
		Bin b = new Bin1DLog10(l, r);
		bs.push_back(b);
	}
	setBins(bs);
}

template<class B>
Histogram1<B>::~Histogram1() {
}

template<class B>
void Histogram1<B>::setBins(vector<Bin> b) {
	bins = b;
	nBins = bins.size();
	contents.resize(nBins, 0.);
}


// Histogram1DLin::Histogram1DLin(double vMin, double vMax, unsigned int n) : Histogram1() {
// 	nBins = n;

// 	vector<Bin> bs;
// 	double dx = (vMax - vMin) / n;
// 	for (size_t i = 0; i < n; i++) {
// 		double l = vMin + i * dx;
// 		double r = vMin + (i + 1) * dx;
// 		Bin b = new Bin1DLin(l, r);
// 		bs.push_back(b);
// 	}
// 	setBins(bs);
// }

// double Histogram1DLin::interpolateAt(const double& x0) const {
// 	if (! isInRange(x0))
// 		throw std::runtime_error("Value is out of range and cannot be interpolated.");

// 	size_t idx = getBinIndex(x0);

// 	vector<double> vx;
// 	vector<double> vy;
// 	for (size_t i = idx - 1; i < idx + 2; i++) {
// 		vx.push_back(bins[i]->getCentre());
// 		vy.push_back(contents[i]);
// 	}

// 	return interpolate(x0, vx, vy);
// }

// std::ostream& operator<<(std::ostream& os, const Histogram1DLin& h) {
// 	os << "Histogram1DLinD: " << endl;
// 	os << "  . number of bins: " << h.getNumberOfBins() << endl;
// 	os << "  . bin edges: " << h.leftEdge() << ", " << h.rightEdge() << endl;
// 	return os;
// }


// Histogram1DLog10::Histogram1DLog10(double vMin, double vMax, unsigned int n) : Histogram1() {
// 	nBins = n;

// 	vector<Bin> bs;
// 	double dx = (log10(vMax) - log10(vMin)) / n;
// 	for (size_t i = 0; i < n; i++) {
// 		double l = pow(10, log10(vMin) + i * dx);
// 		double r = pow(10, log10(vMin) + (i + 1) * dx);
// 		Bin b = new Bin1DLog10(l, r);
// 		bs.push_back(b);
// 	}

// 	setBins(bs);
// }

// double Histogram1DLog10::interpolateAt(const double& x0) const {
// 	if (! isInRange(x0))
// 		throw std::runtime_error("Value is out of range and cannot be interpolated.");

// 	size_t idx = getBinIndex(x0);
// 	vector<double> vx;
// 	vector<double> vy;
// 	for (size_t i = idx - 1; i < idx + 2; i++) {
// 		vx.push_back(log10(bins[i]->getCentre()));
// 		vy.push_back(contents[i]);
// 	}

// 	return interpolate(log10(x0), vx, vy);
// }

// std::ostream& operator<<(std::ostream& os, const Histogram1DLog10& h) {
// 	os << "Histogram1DLog10: " << endl;
// 	os << "  . number of bins: " << h.getNumberOfBins() << endl;
// 	os << "  . bin edges: " << h.leftEdge() << ", " << h.rightEdge() << endl;
// 	return os;
// }








//////////////////

Histogram1D::Histogram1D(string scale, bool logCDF, double minLogCDF) {
	setScale(scale);
	setLogCDF(logCDF);
	setMinimumLogCDF(minLogCDF);
}

Histogram1D::Histogram1D(double vmin, double vmax, unsigned int n, string scale, bool logCDF, double minLogCDF) {
	setScale(scale);
	setLogCDF(logCDF);
	setMinimumLogCDF(minLogCDF);
	setNumberOfBins(n);
	initBins(vmin, vmax, n);
}

Histogram1D::Histogram1D(const Histogram1D& h) {
	reset();

	nBins = h.getNumberOfBins();

	for (size_t i = 0; i < nBins; i++) {
		edges.push_back(h.getBinEdges()[i]);
		centres.push_back(h.getBinCentres()[i]);
		contents.push_back(h.getBinContents()[i]);
	}
	edges.push_back(h.getBinEdges()[nBins]);

	edges.resize(nBins + 1);
	centres.resize(nBins);
	contents.resize(nBins);

	setScale(h.getScale());
	setLogCDF(h.isCDFLogarithmic());
	setMinimumLogCDF(h.getMinimumLogCDF());
}

Histogram1D::Histogram1D(Histogram1D&& h) noexcept {
	*this = std::move(h);
}

Histogram1D::~Histogram1D() {
	// reset();
}

void Histogram1D::initBins(double vmin, double vmax, unsigned int n) {
	if (n <= 0)
		throw std::runtime_error("Number of bins must be greater than zero.");

	reset();

	if (scale == "log10") {
		edges = computeBinEdgesLog10(vmin, vmax, n);
		centres = computeBinCentresLog10(edges, scale); 
	} else {
		edges = computeBinEdgesLinear(vmin, vmax, n);
		centres = computeBinCentresLinear(edges, scale);
	}

	edges.resize(n + 1);
	centres.resize(n);
	contents.resize(n, 0.);
}

void Histogram1D::prepareCDF() {
	if (edges.size() < 2)
		throw std::runtime_error("Number of bin edges must be greater than one.");

	// ensure vector is empty
	cdf.clear();

	// first entry of CDF is 0
	cdf.push_back(0);

	if (logCDF) {
		if (scale == "lin") {
			double cum = 0;
			for (size_t i = 0; i < nBins - 1; i++) {
				cum += contents[i + 1];
				cdf.push_back(cum);
			}
		} else if (scale == "log10") {
			double cum = 0;
			for (size_t i = 1; i < nBins; i++) {
				double dx = edges[i] - edges[i - 1];
				double y = log10(contents[i] + pow(10, minLogCDF - 1));
				// double dy = contents[i] - contents[i - 1];
				// double x = centres[i - 1];
				// double y = contents[i - 1];
				// y = log10(y);
				// double w = y;
				// cum += (y * dx / 2. * w);
				cum += y * dx;
				cdf.push_back(cum);
			}
		}

		double tot = sum();
		for (size_t i = 0; i < nBins; i++) {
			cdf[i] /= tot;
		}

		// // do cdf->log(cdf)
		// cdf[0] = minLogCDF;
		// for (size_t i = 1; i < nBins; i++) {
		// 	cdf[i] = log10(cdf[i]);
		// }



	} else {
		if (scale == "lin") {
			double cum = 0;
			for (size_t i = 0; i < nBins - 1; i++) {
				cum += contents[i + 1];
				cdf.push_back(cum);
			}
		} else if (scale == "log10") {
			double cum = 0;
			for (size_t i = 1; i < nBins; i++) {
				double dx = edges[i] - edges[i - 1];
				double dy = contents[i] - contents[i - 1];
				cum += (contents[i - 1] * dx / 2.);
				cdf.push_back(cum);
			}
		}

		double tot = sum();
		for (size_t i = 0; i < nBins; i++) {
			cdf[i] /= tot;
		}
	}

}	

void Histogram1D::setScale(string s) {
	if (s == "log")
		scale = "log10";
	else
		scale = s;

	if (s != "log" && s != "log10" && s != "lin")
		throw std::runtime_error("Unknown scale " + s + ".");;
}

void Histogram1D::setLogCDF(bool log) {
	logCDF = log;
}

void Histogram1D::setMinimumLogCDF(double v) {
	minLogCDF = v;
}

void Histogram1D::setNumberOfBins(unsigned int n) {
	nBins = n;
}

void Histogram1D::setBinContents(vector<double> values) {
	if (values.size() != contents.size())
		throw runtime_error("Cannot set bin contents: number of values does not match number of bins.");
	contents = values;
}

string Histogram1D::getScale() const {
	return scale;
}

bool Histogram1D::isCDFLogarithmic() const {
	return logCDF;
}

double Histogram1D::getMinimumLogCDF() const {
	return minLogCDF;
}

vector<double> Histogram1D::getBinEdges() const {
	return edges;
}

vector<double> Histogram1D::getBinWidths() const {
	vector<double> binWidths;
	for (size_t i = 0; i < edges.size() - 1; i++)
		binWidths.push_back(edges[i + 1] - edges[i]);
	return binWidths;
}

vector<double> Histogram1D::getBinCentres() const {
	return centres;
}

vector<double> Histogram1D::getBinContents() const {
	return contents;
}

vector<double> Histogram1D::getCDF() const {
	return cdf;
}

size_t Histogram1D::getBinIndex(const double& v) const {
	if  (v < edges.front() or v > edges.back())
		return -1;

	auto bin = std::lower_bound(edges.begin(), edges.end(), v);
	
	// returns the lower edge of the bin (hence -1)
	return bin - edges.begin() - 1;
}

unsigned int Histogram1D::getNumberOfBins() const {
	return nBins;
}

std::pair<double, double> Histogram1D::getEdgesOfBin(const size_t& i) const {
	return std::make_pair(edges[i], edges[i + 1]);
}

double Histogram1D::getBinCentre(const size_t& i) const {
	if (scale == "log10") {
		return pow(10, (log10(edges[i]) + log10(edges[i + 1])) / 2.);
	} else {
		return (edges[i] + edges[i + 1]) / 2.;
	}
}

double Histogram1D::getBinContent(const size_t& i) const {
	return contents[i];
}

void Histogram1D::setBinContent(size_t i, double v) {
	contents[i] = v;
}

void Histogram1D::push(double v, double w) {
	size_t idx = getBinIndex(v);
	if (idx >= 0)
		contents[idx] += w; 			
}

void Histogram1D::fill(const vector<double>& v, const vector<double>& w) {
	if ((w.size() != 0) and (w.size() != v.size()))
		throw runtime_error("Cannot fill histogram: number of weights does not match number of values.");

	for (size_t i = 0; i < v.size(); i++) {
		push(v[i], w.size() == 0 ? 1. : w[i]);
	}
}

void Histogram1D::normalise(double norm) {
	for (size_t i = 0; i < getNumberOfBins(); i++) {
		contents[i] /= norm;
	}
}

double Histogram1D::sum() const {
	return std::accumulate(contents.begin(), contents.end(), 0.);
}

double Histogram1D::integrate() const {
	vector<double> widths = getBinWidths();
	return getIntegral(contents, widths);
}

double Histogram1D::interpolateAt(const double& v) const {
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

double Histogram1D::getSample(bool binned) const {
	// if (cdf.empty())
	// 	throw runtime_error("Cannot get sample: CDF is empty. Call prepareCDF() first.");

	Random& random = Random::instance();

	// return samplingImportance(random);

	// return samplingInverse(random);

	// if (logCDF) {
	// 	return samplingImportance();
	// } else {
	// 	return interpolate(random.rand(), cdf, edges);
	// }

	// if (logCDF) {
	// 	// double r = random.randUniform(minLogCDF, 0.); 
	// 	// r = pow(10, r);
	// 	// double r = random.randExponential();
	// 	double r = random.rand();
	// 	double x = interpolate(r, cdf, edges);
	// 	// size_t idx = getBinIndex(x);
	// 	// return x / centres[idx];
	// 	return x;

	// } else {
	// 	return interpolate(random.rand(), cdf, edges);
	// }
}

double Histogram1D::getSample(const double& xMin, const double& xMax, bool binned) const {
	Random& random = crpropa::Random::instance();

	vector<double> x;
	vector<double> y;

	if (xMin < centres.front() or xMax > centres.back())
		throw runtime_error("Cannot get sample via inversion. Range is out of bounds (considering the bin centres).");

	x.push_back(xMin);
	y.push_back(interpolate(xMin, centres, cdf));
	for (size_t i = 0; i < getNumberOfBins(); i++) {
		if (centres[i] > xMin and centres[i] < xMax) {
			x.push_back(centres[i]);
			y.push_back(cdf[i]);
		}
	}
	x.push_back(xMax);
	y.push_back(interpolate(xMax, centres, cdf));

	return interpolate(random.rand(), y, x);
}

double Histogram1D::getSample(const std::pair<double, double>& range, bool binned) const {
	return getSample(range.first, range.second, binned);
}

double Histogram1D::leftEdge() const {
	return edges.front();
}

double Histogram1D::rightEdge() const {
	return edges.back();
}

void Histogram1D::clear() {
	std::fill(contents.begin(), contents.end(), 0.);
}

void Histogram1D::reset() {
	contents.clear();
	edges.clear();
	centres.clear();
	cdf.clear();
}

double Histogram1D::operator[](const size_t& i) const {
	return getBinContent(i);
}

Histogram1D& Histogram1D::operator=(const Histogram1D& h) {
	if (this == &h)
		return *this;

	reset();
	
	unsigned int n = h.getNumberOfBins();
	for (size_t i = 0; i < n; i++) {
		edges.push_back(h.getBinEdges()[i]);
		centres.push_back(h.getBinCentres()[i]);
		contents.push_back(h.getBinContents()[i]);
	}
	edges.push_back(h.getBinEdges()[n]);

	edges.resize(n + 1);
	centres.resize(n);
	contents.resize(n);

	return *this;
}

vector<double> Histogram1D::computeBinCentres(const vector<double>& binEdges, const string& scale) {
	if (scale == "log10") {
		return computeBinCentresLog10(binEdges, scale);
	} else {
		return computeBinCentresLinear(binEdges, scale);
	}
}

vector<double> Histogram1D::computeBinEdgesLinear(double vmin, double vmax, unsigned int n) {
	vector<double> bins;
	for (size_t i = 0; i <= n; i++) {
		double v = vmin + i * (vmax - vmin) / n;
		bins.push_back(v);
	}
	return bins;
}

vector<double> Histogram1D::computeBinEdgesLog10(double vmin, double vmax, unsigned int n) {
	vector<double> bins = computeBinEdgesLinear(log10(vmin), log10(vmax), n);
	for (size_t i = 0; i <= n; i++) {
		bins[i] = pow(10, bins[i]);
	}

	return bins;
}

vector<double> Histogram1D::computeBinCentresLinear(const vector<double>& binEdges, const string& scale) {
	vector<double> binCentres;
	for (size_t i = 0; i < binEdges.size() - 1; i++) {
		binCentres.push_back((binEdges[i] + binEdges[i + 1]) / 2.);
	}
	return binCentres;
}

vector<double> Histogram1D::computeBinCentresLog10(const vector<double>& binEdges, const string& scale) {
	vector<double> binCentres;
	for (size_t i = 0; i < binEdges.size() - 1; i++) {
		binCentres.push_back(pow(10, (log10(binEdges[i]) + log10(binEdges[i + 1])) / 2.));
	}
	return binCentres;
}

double Histogram1D::getIntegral(const vector<double>& y, const vector<double>& dx) {
	if (y.size() != dx.size())
		throw runtime_error("Cannot compute integral: number of y-values does not match number of bin widths.");

	double integral = 0;
	for (size_t i = 0; i < y.size(); i++) {
		integral += y[i] * dx[i];
	}

	return integral;
}

// double Histogram1D::samplingInverse(Random& random) const {
// 	return interpolate(random.rand(), cdf, edges);
// }

// double Histogram1D::samplingImportance(Random& random, std::function<double(double)> weightFunction) const {

//     // Sample from the proposal distribution
//     double r = random.rand();
//     auto it = std::lower_bound(proposal_cdf.begin(), proposal_cdf.end(), r);
//     size_t idx = std::distance(proposal_cdf.begin(), it) - 1;

//     if (idx >= nBins)
//         idx = nBins - 1;

//     double x1 = edges[idx];
//     double x2 = edges[idx + 1];
//     double y1 = contents[idx];
//     double y2 = contents[idx + 1];

//     if (scale == "log10") {
//         y1 = log10(y1 + 1e-12);
//         y2 = log10(y2 + 1e-12);
//     }


//     double t = (r - proposal_cdf[idx]) / (proposal_cdf[idx + 1] - proposal_cdf[idx]);
//     double x = x1 + t * (x2 - x1);

//     // Compute the importance weight using the provided weight function
//     double proposal_density = (y1 + t * (y2 - y1)) / total;
//     double target_density = weightFunction(x);

//     double weight = target_density / proposal_density;

//     // Accept or reject the sample based on the weight
//     double acceptance = random.rand();
//     if (acceptance <= weight) {
//         return x;
//     } else {
//         // If rejected, recursively call the function to get a new sample
//         return samplingImportance(weightFunction);
//     }
// }


std::ostream& operator<<(std::ostream& os, const Histogram1D& h) {
	os << "Histogram1D: " << endl;
	os << "  . scale: " << h.getScale() << endl;
	os << "  . number of bins: " << h.getNumberOfBins() << endl;
	os << "  . bin edges: " << h.leftEdge() << ", " << h.rightEdge() << endl;
	return os;
}



} // namespace liveprop
