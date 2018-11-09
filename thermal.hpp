
#ifndef _THERMAL_HPP_
#define _THERMAL_HPP_

#include <cmath>

double 	Expansivity(double 	T, double P)
{
	return 	3e-5;
}

double	Conductivity_Kc(double    K298, double 	  T, double 	P)
{
	P = P/1e9;
	
	double 	Klat, Krad;

	// Hofmeister et al. (1999)
	double 	a = 0.3;
	double 	gamma = 1.28;
	double 	alpha = Expansivity(T, P*1e9);
	double 	Kt = 128.1;
	double 	Ktt = 4.6;
	
	Klat = K298*pow(298/T, a)*exp( -(4*gamma+1.0/3.0)*alpha*(T-298) )*( 1 + Ktt*P/Kt );
	Krad = 0.01753 - 0.00010365*T + 2.2451*pow(T,2)/1e7 - 3.407*pow(T,3)/1e11;

	return (Klat + Krad);
}

double 	Conductivity_Kv(double 	T, double Kc, double Ra, double Pr, double lambda)
{
	double 	Kv;
	
	if( Ra>=1e19 )
		Kv = 0.089*Kc*pow(Ra,1.0/3.0);
	else
		Kv = 0.22*Kc*pow(Ra,2.0/7.0)*pow(Pr,-1.0/7.0)*pow(lambda,-3.0/7.0);

	return Kv;
}


#endif
