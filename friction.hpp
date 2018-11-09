
#ifndef _FRICTION_HPP_
#define _FRICTION_HPP_

#include <cmath>

// The function for modeling the power of core-mantle friction
// Only applicable to a core radius of ~350 km
// current_time: reckoned from 4.6 Ga, in Myr

std::vector<double>   InstantFrictionalPower(double   current_time, std::string type)
{
	current_time = current_time/1e3;

	double   a;          		 // Semi-major axis length
	double   aa;         		 // Coefficient
	double   Ie;         		 // Inclination angle
	double   Psigma;     		 // Frictional power

	double   Re = 6371;  		 // Earth radius in km
	double   Psigma0 = 3e20;

    std::vector<double>     results;    results.resize(4);

	if(current_time<=1e-2)
		current_time = 1e-2;
    
	if(type=="nominal")
		a = Re * 41.96 * pow(current_time, 0.2111);
    else if(type=="model2") 	           // Walker et al. (1983)
    {
        double  t = 4.6 - current_time;

    	if(t<=0.6)
    		a = 60.2 * Re * pow( 1 - t/1.4, 2.0/13.0 );
    	else
    		a = 60.2 * Re * pow( 1 - (0.2*t+0.48)/1.4, 2.0/13.0 );
    }
    else if(type=="model3")
    {
        double t = 4.6 - current_time;

        if(t<=0.6)
            a = 60.2 * Re * pow( 1 - t/1.4, 2.0/13.0);
        else if( (t>0.6)&&(t<=3.8) )
            a = 60.2 * Re * pow( 1 - (0.24*t+0.456)/1.4, 2.0/13.0 );
        else
            a = 60.2 * Re * pow( 1 - (0.03*t+1.254)/1.4, 2.0/13.0 );
    }
    else if(type=="model4")     // Model 4 still has some problems
    {
        double  t = 4.6 - current_time;
        a = 60.2 * Re * pow(1 - (1.46/1.4)*pow(1-exp(-t/1.46), -1.0), 2.0/13.0 );
    }	
    
    // Inclination angle
	if( a <= (34.2*Re) )
		Ie = 71.84;
	else
	{
		aa = (a/Re-46.6308)/7.7288;
		Ie =  0.1075*pow(aa,10) - 0.0332*pow(aa,9) - 1.0008*pow(aa,8) + 0.6110*pow(aa,7) + 2.7016*pow(aa,6)
                - 1.7281*pow(aa,5) - 2.3280*pow(aa,4) - 1.4509*pow(aa,3) + 6.9951*pow(aa,2)
                - 6.6208*aa + 5.5828;
    }

    Psigma = Psigma0 * pow(sin(Ie*M_PI/180), 3.0) / pow(a/Re, 9.0/2.0);

    results[0] = current_time;
    results[1] = a;
    results[2] = Ie;
    results[3] = Psigma;
    
    return results;
}


#endif