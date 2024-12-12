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


///////////////////////////////////////////////////////////////////////////////////////////////////

Bin1DLin::Bin1DLin() {
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

bool Bin1DLin::isLinear() const {
	return true;
}

bool Bin1DLin::isLn() const {
	return false;
}

bool Bin1DLin::isLog2() const {
	return false;
}

bool Bin1DLin::isLog10() const {
	return false;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

Bin1DLogarithmic::Bin1DLogarithmic(double b) {
	setBase(b);
}

Bin1DLogarithmic::Bin1DLogarithmic(double l, double r, double b) {
	setEdges(l, r);
	setBase(b);
	setCentre(pow(b, (logBase(l, b) + logBase(r, b)) / 2.));
}

void Bin1DLogarithmic::setBase(double b) {
	base = b;
}

double Bin1DLogarithmic::getBase() const { 
	return base;
}

double Bin1DLogarithmic::rand(Random& random) const {
	double base = getBase();
	double vMin = logBase(getLeftEdge(), base);
	double vMax = logBase(getRightEdge(), base);
	
	double r = random.randUniform(vMin, vMax);

	return pow(base, r);
}

double Bin1DLogarithmic::directTransformation(const double& v) const {
	return logBase(v, getBase());
}

double Bin1DLogarithmic::inverseTransformation(const double& v) const {
	return pow(getBase(), v);
}

bool Bin1DLogarithmic::isLinear() const {
	return false;
}

bool Bin1DLogarithmic::isLn() const {
	return (base == exp(1.));
}

bool Bin1DLogarithmic::isLog2() const {
	return (base == 2.);
}

bool Bin1DLogarithmic::isLog10() const {
	return (base == 10.);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

Bin1DLog10::Bin1DLog10() {
	setBase(10.);
}

Bin1DLog10::Bin1DLog10(double l, double r) : Bin1DLogarithmic(l, r, 10.) {
}

Bin1DLog2::Bin1DLog2() {
	setBase(2.);
}

Bin1DLog2::Bin1DLog2(double l, double r) : Bin1DLogarithmic(l, r, 2.) {
}

Bin1DLn::Bin1DLn() {
	setBase(exp(1.));
}

Bin1DLn::Bin1DLn(double l, double r) : Bin1DLogarithmic(l, r, exp(1.)) {
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

void Histogram1D::setIsPDF(bool b) {
	isPDF = b;
}

void Histogram1D::setIsCDF(bool b) {
	isCDF = b;
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

vector<Histogram1D::Bin> Histogram1D::getBins() const {
	return bins;
}

bool Histogram1D::getIsPDF() const {
	return isPDF;
}

bool Histogram1D::getIsCDF() const {
	return isCDF;
}

bool Histogram1D::isLinear() const {
	return bins[0]->isLinear();
}

bool Histogram1D::isLog10() const {
	return bins[0]->isLog10();
}

bool Histogram1D::isLog2() const {
	return bins[0]->isLog2();
}

bool Histogram1D::isLn() const {
	return bins[0]->isLn();
}

bool Histogram1D::isIrregular() const {
	return ! isRegular();
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

vector<double> Histogram1D::computeVectorCDF() const {
	vector<double> cdf = {0.};

	double cum = 0;
	for (size_t i = 1; i < nBins; i++) {
		auto bin = getBin(i);
		double dx = bin->getWidth();
		double y = contents[i];
		cum += y * dx;
		cdf.push_back(cum);
	}

	for (size_t i = 1; i < nBins; i++) {
		cdf[i] /= cdf.back();
	}

	return cdf;
}

void Histogram1D::makePDF() {
	if (isPDF)
		return;

	double integral = integrate();
	if (integral != 1.)
		normalise(integral);

	isPDF = true;
	isCDF = false;
}

void Histogram1D::makeCDF() {
	if (isCDF)
		return;

	contents = computeVectorCDF();
	isCDF = true;
	isPDF = false;
}

void Histogram1D::reset() {
	std::fill(contents.begin(), contents.end(), 0.);
	std::fill(weights.begin(), weights.end(), 1.);
}

void Histogram1D::clear() {
	contents.clear();
	bins.clear();
	weights.clear();
}

double Histogram1D::operator[](const size_t& i) const {
	return contents[i];
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void RegularHistogram1D::setBins(vector<Bin> b) {
	bins = b;
	nBins = bins.size();
	contents.resize(nBins, 0.);
	weights.resize(nBins, 0.);
}

bool RegularHistogram1D::isRegular() const {
	return true;
}

double RegularHistogram1D::directTransformation(const double& v) const {
	return bins[0]->directTransformation(v);
}

double RegularHistogram1D::inverseTransformation(const double& v) const {
	return bins[0]->inverseTransformation(v);
}

double RegularHistogram1D::interpolateAt(const double& x0) const {
	if (! isInRange(x0)) {
		cout << "Value out of range for interpolation. Returning NaN." << endl;
		return std::nan("");
	}

	std::function<double(double)> interp = getInterpolator();

	return interp(x0);
}

std::function<double(double)> RegularHistogram1D::getInterpolator(const std::pair<double, double>& range) const {
	return [this](const double& x0) { 
		// get index of bin containing the value
		size_t idx = getBinIndex(x0);

		// value is in the first bin, within the histogram range, but below the first bin centre
		if (x0 < bins[0]->getCentre()) {
			vector<double> x;
			vector<double> y;
			for (size_t i = 1; i < 3; i++) {
				x.push_back(directTransformation(bins[i]->getCentre()));
				y.push_back(contents[i]);
			}
			x[0] = directTransformation(bins[0]->getLeftEdge());
			y[0] = twoPointExtrapolation(x[0], x[1], y[1], x[2], y[2]);
			return interpolate(directTransformation(x0), x, y);

		// value is in the last bin, within the histogram range, but above the last bin centre
		} else if (x0 >= bins[nBins - 1]->getCentre()) {
			vector<double> x;
			vector<double> y;
			for (size_t i = 0; i < 2; i++) {
				size_t j = nBins - 3 + i + 1;
				x.push_back(directTransformation(bins[j]->getCentre()));
				y.push_back(contents[j]);
			}
			size_t j = nBins - 2;
			x[2] = directTransformation(bins[j]->getRightEdge());
			y[2] = twoPointExtrapolation(x[2], x[0], y[0], x[1], y[1]);
			return interpolate(directTransformation(x0), x, y);
		}

		return interpolateEquidistant(directTransformation(x0), directTransformation(leftEdge()), directTransformation(rightEdge()), contents);
	};
}

ref_ptr<Histogram1D> RegularHistogram1D::getHistogramPDF() const {
	if (isPDF)
		return const_cast<RegularHistogram1D*>(this);
	
	ref_ptr<Histogram1D> h = clone();
	h->makePDF();

	return h;
}

ref_ptr<Histogram1D> RegularHistogram1D::getHistogramCDF() const {
	if (isCDF)
		return const_cast<RegularHistogram1D*>(this);
	
	ref_ptr<Histogram1D> h = clone();
	h->makeCDF();

	return h;
}

// ref_ptr<Histogram1D>& RegularHistogram1D::operator=(const RegularHistogram1D& h) {
// 	if (this == &h) {
// 		ref_ptr<Histogram1D> hNew = *this;
// 		return hNew;
// 	}

// 	nBins = h.getNumberOfBins();
// 	isPDF = h.isPDF;
// 	isCDF = h.isCDF;

// 	setBins(h.getBins());
// 	setBinContents(h.getBinContents());

// 	return *this;
// }

std::ostream& operator<<(std::ostream& os, const RegularHistogram1D& h) {
	os << "Histogram1D " << endl;
	os << ". number of bins: " << h.getNumberOfBins();

	if (h.isLinear()) {
		os << " (linear spacing)" << endl;
	} else if (h.isLog10()) {
		os << " (log10 spacing)" << endl;
	} else if (h.isLn()) {
		os << " (ln spacing)" << endl;
	} else if (h.isLog2()) {
		os << " (log2 spacing)" << endl;
	}

	os << ". edges: " << h.leftEdge() << ", " << h.rightEdge() << endl;
	vector<double> v = h.getBinContents();

	os << ". extrema: " << *std::max_element(v.begin(), v.end()) << ", " <<  *std::max_element(v.begin(), v.end()) << endl;

	return os;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

Histogram1DLin::Histogram1DLin() {
	nBins = 0;
	isPDF = false;
	isCDF = false;
}

Histogram1DLin::Histogram1DLin(double vMin, double vMax, unsigned int n) {
	nBins = n;
	isPDF = false;
	isCDF = false;

	vector<Bin> bins;
	double dx = (vMax - vMin) / n;
	for (size_t i = 0; i < n; i++) {
		double l = vMin + i * dx;
		double r = l + dx;
		Bin bin = new Bin1DLin(l, r);
		bins.push_back(bin);
	}

	setBins(bins);
}

ref_ptr<Histogram1D> Histogram1DLin::clone() const {
	ref_ptr<Histogram1D> h = new Histogram1DLin(leftEdge(), rightEdge(), getNumberOfBins());
	return h;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

Histogram1DLog10::Histogram1DLog10() {
	nBins = 0;
	isPDF = false;
	isCDF = false;
}

Histogram1DLog10::Histogram1DLog10(double vMin, double vMax, unsigned int n) {
	nBins = n;
	isPDF = false;
	isCDF = false;

	vector<Bin> bins;
	Bin bin = new Bin1DLog10();
	double dx = (bin->directTransformation(vMax) - bin->directTransformation(vMin)) / n;

	for (size_t i = 0; i < n; i++) {
		double l = bin->directTransformation(vMin) + i * dx;
		double r = l + dx;
		l = bin->inverseTransformation(l);
		r = bin->inverseTransformation(r);
		bin = new Bin1DLog10(l, r);
		bins.push_back(bin);
	}

	setBins(bins);
}

ref_ptr<Histogram1D> Histogram1DLog10::clone() const {
	ref_ptr<Histogram1D> h = new Histogram1DLog10(leftEdge(), rightEdge(), getNumberOfBins());
	return h;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

Histogram1DLog2::Histogram1DLog2() {
	nBins = 0;
	isPDF = false;
	isCDF = false;
}

Histogram1DLog2::Histogram1DLog2(double vMin, double vMax, unsigned int n) {
	nBins = n;
	isPDF = false;
	isCDF = false;
	vector<Bin> bins;
	Bin bin = new Bin1DLog2(vMin, vMax);
	double dx = (bin->directTransformation(vMax) - bin->directTransformation(vMin)) / n;

	for (size_t i = 0; i < n; i++) {
		double l = bin->directTransformation(vMin) + i * dx;
		double r = l + dx;
		l = bin->inverseTransformation(l);
		r = bin->inverseTransformation(r);
		bin = new Bin1DLog2(l, r);
		bins.push_back(bin);
	}

	setBins(bins);
}

ref_ptr<Histogram1D> Histogram1DLog2::clone() const {
	ref_ptr<Histogram1D> h = new Histogram1DLog2(leftEdge(), rightEdge(), getNumberOfBins());
	return h;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

Histogram1DLn::Histogram1DLn() {
	nBins = 0;
	isPDF = false;
	isCDF = false;
}

Histogram1DLn::Histogram1DLn(double vMin, double vMax, unsigned int n) {
	nBins = n;
	isPDF = false;
	isCDF = false;

	vector<Bin> bins;
	Bin bin = new Bin1DLn(vMin, vMax);
	double dx = (bin->directTransformation(vMax) - bin->directTransformation(vMin)) / n;

	for (size_t i = 0; i < n; i++) {
		double l = bin->directTransformation(vMin) + i * dx;
		double r = l + dx;
		l = bin->inverseTransformation(l);
		r = bin->inverseTransformation(r);
		bin = new Bin1DLog10(l, r);
		bins.push_back(bin);
	}

	setBins(bins);
}

ref_ptr<Histogram1D> Histogram1DLn::clone() const {
	ref_ptr<Histogram1D> h = new Histogram1DLn(leftEdge(), rightEdge(), getNumberOfBins());
	return h;
}



} // namespace livpropa
