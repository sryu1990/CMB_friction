
/* A header file to model the solidus, liquidus and rheoligy transition in solid mantle.


Shuoran Yu

State Key Laboratory of Lunar and Planetary Sciences
Macau Unviersity of Science and Technology

E-mail: shuoran.yu@icloud.com
GitHub homepage: https://github.com/sryu1990

*/



#ifndef _MELT_HPP_
#define _MELT_HPP_

#include <cmath>

class _Melt
{
	private:

		double 	Tsol0, Tsol1, Tsol2, Tsol3;
		double 	Tliq0, Tliq1, Tliq2, Tliq3;
		double 	phi_0;

		double 	rhom;
		double 	rhoc;
		double 	rc;
		double 	rp;

		double 	G;

	public:

		_Melt(double rhom_in, double rhoc_in, double rp_in, double rc_in)
		{
			G = 6.67e-11;

			rhom = rhom_in;
			rhoc = rhoc_in;
			rc = rc_in;
			rp = rp_in;
			
			// Threshold melt fraction for rheological transition
			phi_0 = 0.60;
			
			// Data from Laneuville et al. (2013), JGR, for peridotite
			Tsol0 = 1409;  Tsol1 = 134.2;  Tsol2 = -6.581;  Tsol3 = 0.1054;
			Tliq0 = 2035;  Tliq1 = 57.46;  Tliq2 = -3.4872; Tsol3 = 0.0769;
		};

		double  getGravity(double r)
		{
			double gr = 4*M_PI*G*(r*rhom + pow(rc,3)*(rhoc-rhom)/pow(r,2))/3;
			return gr;
		};

		double 	getPressure(double r)
		{
			double pr = 4*M_PI*rhom*G*pow(rc,3)*(rhoc-rhom)*(1/r-1/rp)/3 + 2*M_PI*G*pow(rhom,2)*(pow(rp,2)-pow(r,2))/3;
			return pr;
		};

		double 	getSolidusTemp(double r)
		{
			double p = getPressure(r);
			p = p/1e9;

			double Tsol = Tsol0 + Tsol1 * p + Tsol2 * pow(p,2) + Tsol3 * pow(p,3);
			return Tsol;
		};

		double getSolidusTempGradient(double r)
		{
			double p = getPressure(r);
			p = p/1e9;

			double 	dTsoldp = Tsol1 + 2.0*Tsol2*p + 3.0*Tsol3*pow(p,2);
			double 	dpdr = - 4 * M_PI * pow(r,2) * rhom * getGravity(r);

			return dTsoldp*dpdr;
		};

		double getLiquidusTemp(double r)
		{
			double p = getPressure(r);
			p = p/1e9;

			double Tliq = Tliq0 + Tliq1 * p + Tliq2 * pow(p,2) + Tliq3 * pow(p,3);
			return Tliq;
		};

		double getLiquidusTempGradient(double r)
		{
			double p = getPressure(r);
			p = p/1e9;

			double 	dTliqdp = Tliq1 + 2.0*Tliq2*p + 3.0*Tliq3*pow(p,2);
			double 	dpdr = - 4 * M_PI * pow(r,2) * rhom * getGravity(r);

			return dTliqdp*dpdr;
		};

		double getRheologyTemp(double r)
		{
			double 	Tsol = getSolidusTemp(r);
			double 	Tliq = getLiquidusTemp(r);

			return Tsol*phi_0 + Tliq*(1-phi_0);
		};

		double getRheologyTempGradient(double r)
		{
			double 	dTliqdr = getLiquidusTempGradient(r);
			double 	dTsoldr = getSolidusTempGradient(r);

			return  (1-phi_0)*dTsoldr + phi_0*dTliqdr;
		};

		double getMeltFraction(double Temp, double r)
		{
			// A linear melting model
			double frac = ( Temp - getSolidusTemp(r) )/( getLiquidusTemp(r) - getSolidusTemp(r) );

			if( frac<=0.0 )
				frac = 0.0;
			else if( frac>=1.0 )
				frac = 1.0;

			return frac;
		}

		double getViscosity(double T, double r)
		{
			double 	eta_l, eta_s, eta, phi;
			double 	E, V, eta_r, Tr, a, R;
			double 	A, B, C;

			// Data for peridotite
			E = 375e3;  V = 30e-6;  eta_r = 1e21;  Tr = 1600;  a = 26;  R = 8.314;
			A = 

			phi = getMeltFraction(T, r);
			// eta_l = exp( A + B/(T-C) );		// Viscosity of magma
			eta_l = 1e3;
			eta_s = eta_r * exp( E/(R*T) - E/(R*Tr) );	// Arrhnenius law for diffusion creep

			if(phi>phi_0)	// Einstein-Roscoe formula
				eta = eta_l / pow(phi, 2.5);
			else
				eta = eta_s * exp(-a*phi);

			return eta;
		}
};



#endif