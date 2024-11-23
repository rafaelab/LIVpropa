#include "livpropa/Histogram.h"

namespace livpropa {


///////////////////////////////////////////////////////////////////////////////////////////////////

void Bin1D::setEdges(double l, double r) {
	left = l;
	right = r;
}

void Bin1D::setCentre(double c) {
	centre = c;
}

pair<double, double> Bin1D::getEdges() const {
	return std::make_pair(left, right);
}

double Bin1D::getLeftEdge() const {
	return left;
}

double Bin1D::getRightEdge() const {
	return right;
}

double Bin1D::getCentre() const {
	return centre;
}

double Bin1D::getWidth() const {
	return right - left;
}

bool Bin1D::isInBin(const double& v) const {
	return v >= left and v < right;
}


Bin1DLin::Bin1DLin(double l, double r) {
	setEdges(l, r);
	setCentre((l + r) / 2.);
}

double Bin1DLin::rand(Random& random) const {
	return random.randUniform(getLeftEdge(), getRightEdge());
}

template<LogBase B>
Bin1DLog<B>::Bin1DLog(double l, double r) {
	setEdges(l, r);

	double b = getBase();
	setCentre(pow(b, (logBase(l, b) + logBase(r, b)) / 2.));
}

template<LogBase B>
double Bin1DLog<B>::getBase() const { 
	switch (B) {
		case LogBase::e:
			return M_E;
		case LogBase::two:
			return 2.;
		case LogBase::ten:
			return 10.;
	}
}

template<LogBase B>
double Bin1DLog<B>::rand(Random& random) const {
	return random.randUniform(getLeftEdge(), getRightEdge());
}


///////////////////////////////////////////////////////////////////////////////////////////////////

bool Histogram1D::isInRange(const double& v) const {
	if (v >= bins[0]->getLeftEdge() and v < bins[nBins - 1]->getRightEdge())
		return true;

	return false;
}

size_t Histogram1D::getBinIndex(const double& v) const {
	if (! isInRange(v))
		return -1;

	for (size_t i = 0; i < nBins; i++) {
		if (bins[i]->isInBin(v))
			return i;
	}

	return -1;
}

unsigned int Histogram1D::getNumberOfBins() const {
	return nBins;
}

double Histogram1D::leftEdge() const {
	return bins[0]->getLeftEdge();
}

double Histogram1D::rightEdge() const {
	return bins[nBins - 1]->getRightEdge();
}

Histogram1D::Bin Histogram1D::getBin(const size_t& i) const {
	return bins[i];
}

double Histogram1D::getBinContent(const size_t& i) const {
	return contents[i];
}

double Histogram1D::getBinCentre(const size_t& i) const {
	return bins[i]->getCentre();
}

vector<double> Histogram1D::getBinEdges() const {
	vector<double> edges;
	for (size_t i = 0; i < nBins; i++) {
		edges.push_back(bins[i]->getLeftEdge());
	}
	edges.push_back(bins[nBins - 1]->getRightEdge());

	return edges;
}

vector<double> Histogram1D::getBinCentres() const {
	vector<double> centres;
	for (size_t i = 0; i < nBins; i++) {
		centres.push_back(bins[i]->getCentre());
	}

	return centres;
}

vector<double> Histogram1D::getBinContents() const {
	return contents;
}

void Histogram1D::push(const double& v, const double& w) {
	size_t i = getBinIndex(v);
	contents[i] += w;
	weights[i] = w;
}

void Histogram1D::fill(const vector<double>& v, const vector<double>& w) {
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

void Histogram1D::setBinContent(const size_t& idx, const double& value) {
	contents[idx] = value;
}

void Histogram1D::setBinContents(const vector<double>& values) {
	if (values.size() != nBins)
		throw std::runtime_error("Number of values must match number of bins.");

	contents = values;
}

void Histogram1D::normalise(double norm) {
	for (size_t i = 0; i < nBins; i++) {
		contents[i] /= norm;
	}
}

double Histogram1D::sum() const {
	double s = 0;
	for (size_t i = 0; i < nBins; i++) {
		s += contents[i];
	}

	return s;
}

double Histogram1D::integrate() const {
	double s = 0;
	for (size_t i = 0; i < nBins; i++) {
		double dx = bins[i]->getWidth();
		s += contents[i] * dx;
	}

	return s;
}

double Histogram1D::operator[](const size_t& i) const {
	return contents[i];
}

void Histogram1D::reset() {
	contents.clear();
	bins.clear();
	weights.clear();
}

void Histogram1D::clear() {
	std::fill(contents.begin(), contents.end(), 0.);
	std::fill(weights.begin(), weights.end(), 1.);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

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
Histogram1<B>::Histogram1(const Histogram1<B>& h) {
	reset();

	nBins = h.getNumberOfBins();
	vector<Bin> newBins;
	
	for (size_t i = 0; i < nBins; i++) {
		Bin b0 = h.getBin(i);
		Bin bin = new B(b0->getLeftEdge(), b0->getRightEdge());
		newBins.push_back(bin);
		setBinContent(i, h.getBinContent[i]);
	}
}

template<class B>
Histogram1<B>::Histogram1(Histogram1<B>&& h) noexcept {
	*this = std::move(h);
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

// template<>
// double Histogram1<Bin1DLin>::interpolateAt(const double& v) const {
// 	Random& random = Random::instance();

// 	size_t i = getBinIndex(v);

// 	if (binned) {
// 		size_t idx = random.randBin(v) + 1;
// 		return random.randUniform();
// 	}

// }



} // namespace liveprop
