/*
 * Koordinatentransformation.cpp
 *
 *  Created on: 08.09.2010
 *      Author: martin
 */

#include<math.h>
#include<iostream>


using namespace std;

////////////////////////////////////////////////////////////////////////////////////////
// Kugelkoordinaten
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////
// Funktionsstart int  Umwandlung_Kugel_in_Karthesisch(double r,double phi,double theta,double& x,double& y,double& z);
////////////////////////////////////////////////
void  Umwandlung_Kugel_in_Karthesisch(double r,double phi,double theta,double& x,double& y,double& z)
{
    // Winkel in Grad
    const double pi=3.1415926535897;

    x=r*cos(theta*pi/180.0)*cos(phi*pi/180.0);
    y=r*cos(theta*pi/180.0)*sin(phi*pi/180.0);
    z=r*sin(theta*pi/180.0);
}
////////////////////////////////////////////////
// ENDE int  Umwandlung_Kugel_in_Karthesisch(double r,double phi,double theta,double& x,double& y,double& z);
////////////////////////////////////////////////

////////////////////////////////////////////////
// Funktionsstart int Umwandlung_Karthesisch_in_Kugel(double x,double y,double z,double& r,double& phi,double& theta);
////////////////////////////////////////////////
void Umwandlung_Karthesisch_in_Kugel(double x,double y,double z,double& r,double& phi,double& theta)
{
    // Winkel in Grad
    // Der Punkt 0,0,0 ist ausgeschlossen
    const double pi=3.1415926535897;
    const double epsilon=0.001; //

    r=sqrt(x*x+y*y+z*z);
    if (r==0)
    {    cout<<" Umwandlung in Kugelkoordinaten von (0,0,0) nicht sinnvoll\n";
         phi=0; theta=0;
         return;
    }
    theta=180.0/pi*asin(z/r);
    if ((theta+epsilon>180.0) || (theta-epsilon<-180.0))
    {
        // Phi quasi beliebig
        phi=0;
        return;
    }

    if(((x-epsilon<0)) && ((x+epsilon)>0))
    {
        if (y>0)
            phi=90.0;
        else
            phi=-90.0;
        return;
    }
    phi=180.0/pi*atan2(y,x);
}
////////////////////////////////////////////////
// ENDE int Umwandlung_Karthesisch_in_Kugel(double x,double y,double z,double& r,double& phi,double& theta);
////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
// ENDE  Kugelkoordinaten
////////////////////////////////////////////////////////////////////////////////////////
