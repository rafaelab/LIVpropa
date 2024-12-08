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

double Bin1D::randUniform(Random& random) const {
	return random.randUniform(left, right);
}



Bin1DLin::Bin1DLin(double l, double r) {
	setEdges(l, r);
	setCentre((l + r) / 2.);
}

double Bin1DLin::rand(Random& random) const {
	return random.randUniform(getLeftEdge(), getRightEdge());
}

double Bin1DLin::directTransformation(const double& v) const {
	return v;
}

double Bin1DLin::inverseTransformation(const double& v) const {
	return v;
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
	double base = getBase();
	double vMin = logBase(getLeftEdge(), base);
	double vMax = logBase(getRightEdge(), base);
	
	double r = random.randUniform(vMin, vMax);

	return pow(base, r);
}

template<LogBase B>
double Bin1DLog<B>::directTransformation(const double& v) const {
	return logBase(v, getBase());
}

template<LogBase B>
double Bin1DLog<B>::inverseTransformation(const double& v) const {
	return pow(getBase(), v);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

bool Histogram1D::isInRange(const double& v) const {
	if (v >= bins[0]->getLeftEdge() and v < bins[nBins - 1]->getRightEdge())
		return true;

	return false;
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

void Histogram1D::setBinContent(const size_t& idx, const double& value) {
	contents[idx] = value;
}

void Histogram1D::setBinContents(const vector<double>& values) {
	if (values.size() != nBins)
		throw std::runtime_error("Number of values must match number of bins.");

	contents = values;
}

Histogram1D::Bin Histogram1D::getBin(const size_t& i) const {
	return bins[i];
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

double Histogram1D::getBinContent(const size_t& i) const {
	return contents[i];
}

double Histogram1D::getBinCentre(const size_t& i) const {
	return bins[i]->getCentre();
}

double Histogram1D::getBinWidth(const size_t& i) const {
	return bins[i]->getWidth();
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

template<class B>
double Histogram1<B>::directTransformation(const double& v) const {
	return bins[0]->directTransformation(v);
}

template<class B>
double Histogram1<B>::inverseTransformation(const double& v) const {
	return bins[0]->inverseTransformation(v);
}

template<class B>
double Histogram1<B>::interpolateAt(const double& x0) const {
	if (! isInRange(x0)) {
		cout << "Value out of range for interpolation. Returning NaN." << endl;
		return std::nan("");
	}

	vector<double> x;
	vector<double> y;

	size_t idx = getBinIndex(x0);

	// value is in the first bin, within the histogram range, but below the first bin centre
	if (x0 < bins[0]->getCentre()) {
		for (size_t i = 1; i < 3; i++) {
			x.push_back(directTransformation(bins[i]->getCentre()));
			y.push_back(contents[i]);
		}
		x[0] = directTransformation(bins[0]->getLeftEdge());
		y[0] = twoPointExtrapolation(x[0], x[1], y[1], x[2], y[2]);

	// value is in the last bin, within the histogram range, but above the last bin centre
	} else if (x0 >= bins[nBins - 1]->getCentre()) {
		for (size_t i = 0; i < 2; i++) {
			size_t j = nBins - 3 + i + 1;
			x.push_back(directTransformation(bins[j]->getCentre()));
			y.push_back(contents[j]);
		}
		size_t j = nBins - 2;
		x[2] = directTransformation(bins[j]->getRightEdge());
		y[2] = twoPointExtrapolation(x[2], x[0], y[0], x[1], y[1]);

	// general case
	} else {
		for (size_t i = idx - 1; i < idx + 2; i++) {
			x.push_back(directTransformation(bins[i]->getCentre()));
			y.push_back(contents[i]);
		}
	}

	return interpolate(directTransformation(x0), x, y);
}

template<class B>
std::ostream& operator<<(std::ostream& os, const Histogram1<B>& h) {
	os << "Histogram1D " << endl;
	os << ". number of bins: " << h.getNumberOfBins();

	if (std::is_same<B, Bin1DLin>::value) {
		os << " (linear spacing)" << endl;
	} else if (std::is_same<B, Bin1DLog10>::value) {
		os << " (log10 spacing)" << endl;
	}

	os << ". edges: " << h.leftEdge() << ", " << h.rightEdge() << endl;
	vector<double> v = h.getBinContents();

	os << ". extrema: " << *std::max_element(v.begin(), v.end()) << ", " <<  *std::max_element(v.begin(), v.end()) << endl;

	return os;
}


} // namespace livpropa
