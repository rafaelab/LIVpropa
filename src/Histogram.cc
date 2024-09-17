#include "livpropa/Histogram.h"

namespace livpropa {

Histogram1D::Histogram1D(string scale) {
	setScale(scale);
}

Histogram1D::Histogram1D(double vmin, double vmax, unsigned int n, string scale) {
	setScale(scale);
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

void Histogram1D::setScale(string s) {
	if (s == "log")
		scale = "log10";
	else
		scale = s;

	if (s != "log" && s != "log10" && s != "lin")
		throw std::runtime_error("Unknown scale " + s + ".");;
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

double Histogram1D::getSample() const {
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

double Histogram1D::getSample(const double& xMin, const double& xMax) const {
	double vMin = (xMin < edges.front()) ? edges.front() : xMin;
	double vMax = (xMax > edges.back()) ? edges.back() : xMax;

	size_t idxL = getBinIndex(vMin);
	size_t idxU = getBinIndex(vMax) + 1;

	// possible bugs in first/last bins!

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


std::ostream& operator<<(std::ostream& os, const Histogram1D& h) {
	os << "Histogram1D: " << endl;
	os << "  . scale: " << h.getScale() << endl;
	os << "  . number of bins: " << h.getNumberOfBins() << endl;
	os << "  . bin edges: " << h.getBinEdges().front() << ", " << h.getBinEdges().back() << endl;
	return os;
}





} // namespace livpropa
