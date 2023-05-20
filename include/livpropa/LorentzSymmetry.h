#ifndef LIVPROPA_LORENTZSYMMETRY_H
#define LIVPROPA_LORENTZSYMMETRY_H

#include <cstdlib>
#include <cstring>
#include <fstream>


#include <crpropa/Units.h>



namespace livpropa {


class LorentzSymmetry {
	protected:
		unsigned int order;
		double parameter;
		double energyScale;
	public:
		LorentzSymmetry(double energy, unsigned int order, double parameter);
		void setOrder(unsigned int order);
		void setParameter(double parameter);
		void setEnergyScale(double energy);
		unsigned int getOrder() const;
		double getParameter() const;
		double getEnergyScale() const;
		bool isSuperluminal() const;
		bool isSubluminal() const;

};



} // namespace livpropa

#endif // LIVPROPA_LORENTZSYMMETRY_H
