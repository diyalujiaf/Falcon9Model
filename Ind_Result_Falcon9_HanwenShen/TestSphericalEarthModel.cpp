#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SphericalEarthModel.h"
#include "SphericalEarthModelConstants.h"
#include "SphericalEarthModel.cpp"

const double DEGTORAD = atan(1.0) / 45.0;
const double RADTODEG = 1.0 / DEGTORAD;

int main(int argc, char **argv)
{
   SphericalEarthModel earth;

   //   Approximate lat/lon/alt of classroom.

   double latitude  =  39.948775 * DEGTORAD;
   double longitude = -75.123262 * DEGTORAD;
   double altitude  =   5.0;

   double xECEF = 0.0;
   double yECEF = 0.0;
   double zECEF = 0.0;

   earth.ComputeECEFFromLLA(latitude, longitude, altitude, xECEF, yECEF, zECEF);

   double xECEFM =  1.251152641869611e+06;
   double yECEFM = -4.709871851986261e+06;
   double zECEFM =  4.081698390505012e+06;

   double deltaX = xECEF - xECEFM;
   double deltaY = yECEF - yECEFM;
   double deltaZ = zECEF - zECEFM;

   printf("Testing ComputeECEFFromLLA: \n");
   printf("Lat: %8.3f, Lon: %8.3f, Alt: %8.3f\n", latitude * RADTODEG, longitude * RADTODEG, altitude);
   printf("xECEF: %8.3f, yECEF: %8.3f, zECEF: %8.3f\n", xECEF, yECEF, zECEF);
   printf("Difference between MATLAB and C++:\n");
   printf("deltaX: %30.20e, deltaY: %30.20e, deltaZ: %30.20e\n", deltaX, deltaY, deltaZ);

   double latitude2  = 0.0;
   double longitude2 = 0.0;
   double altitude2  = 0.0;

   earth.ComputeLLAFromECEF(xECEF, yECEF, zECEF, latitude2, longitude2, altitude2);

   double deltaLat = latitude2 - latitude;
   double deltaLon = longitude2 - longitude;
   double deltaAlt = altitude2 - altitude;

   printf("Testing ComputeLLAFromECEF: \n");
   printf("Lat2: %8.3f, Lon2: %8.3f, Alt2: %8.3f\n", latitude2 * RADTODEG, longitude2 * RADTODEG, altitude2);
   printf("deltaLat: %12.8f, deltaLon: %12.8f, deltaAlt: %12.8f\n", deltaLat * RADTODEG, deltaLon * RADTODEG, deltaAlt);

   double enuLatitude  =  39.948775 * DEGTORAD;
   double enuLongitude = -75.123262 * DEGTORAD;
   double enuAltitude  =   5.0;

   double xENU = 0.0;
   double yENU = 0.0;
   double zENU = 0.0;

   bool positionData = true;

   double xECEF2 = 0.0;
   double yECEF2 = 0.0;
   double zECEF2 = 0.0;

   earth.ConvertENUToECEF(enuLatitude, enuLongitude, enuAltitude, xENU, yENU, zENU, positionData, xECEF2, yECEF2, zECEF2);

   deltaX = xECEF2 - xECEF;
   deltaY = yECEF2 - yECEF;
   deltaZ = zECEF2 - zECEF;

   printf("Testing ConvertENUToECEF: \n");
   printf("xECEF2: %8.3f, yECEF2: %8.3f, zECEF2: %8.3f\n", xECEF2, yECEF2, zECEF2);
   printf("deltaX: %30.20e, deltaY: %30.20e, deltaZ: %30.20e\n", deltaX, deltaY, deltaZ);

   double gX = 0.0;
   double gY = 0.0;
   double gZ = 0.0;

   earth.ComputeGravitationalAcceleration(xECEF, yECEF, zECEF, gX, gY, gZ);

   double gXM = -1.930161831179626;
   double gYM =  7.265951870482143;
   double gZM = -6.296864328214498;

   double deltaGX = gX - gXM;
   double deltaGY = gY - gYM;
   double deltaGZ = gZ - gZM;

   printf("Testing ComputeGravitationalAcceleration: \n");
   printf("gX: %8.3f, gY: %8.3f, gZ: %8.3f\n", gX, gY, gZ);
   printf("Difference between MATLAB and C++:\n");
   printf("deltaGX: %30.20e, deltaGY: %30.20e, deltaGZ: %30.20e\n", deltaGX, deltaGY, deltaGZ);

   double airDensity       = 0.0;
   double speedOfSound     = 0.0;
   double pressure         = 0.0;
   double temperature      = 0.0;
   double dynamicViscosity = 0.0;

   earth.ComputeAtmosphere(altitude, airDensity, speedOfSound, pressure, temperature, dynamicViscosity);

   double airDensityM       = 1.224411247820207;
   double speedOfSoundM     = 3.402808085143658e+02;
   double pressureM         = 1.012649487717776e+05;
   double temperatureM      = 2.881175000255633e+02;
   double dynamicViscosityM = 1.789223457824902e-05;

   double deltaRho = airDensity - airDensityM;
   double deltaSOS = speedOfSound - speedOfSoundM;
   double deltaP   = pressure - pressureM;
   double deltaT   = temperature - temperatureM;
   double deltaDV  = dynamicViscosity - dynamicViscosityM;

   printf("Testing ComputeAtmosphere: \n");
   printf("airDensity: %8.3f, speedOfSound: %8.3f, pressure: %8.3f, temperature: %8.3f\n",
           airDensity, speedOfSound, pressure, temperature);
   printf("Difference between MATLAB and C++:\n");
   printf("deltaRho: %30.20e\n", deltaRho);
   printf("deltaSOS: %30.20e\n", deltaSOS);
   printf("deltaP: %30.20e\n", deltaP);
   printf("deltaT: %30.20e\n", deltaT);
   printf("deltaDV: %30.20e\n", deltaDV);

   double nHatX = 0.0;
   double nHatY = 0.0;
   double nHatZ = 0.0;

   earth.ComputeOutwardUnitNormal(latitude, longitude, nHatX, nHatY, nHatZ);

   double nHatXM =  0.196822040918198;
   double nHatYM = -0.740922058067887;
   double nHatZM =  0.642102474747794;

   double deltaNX = nHatX - nHatXM;
   double deltaNY = nHatY - nHatYM;
   double deltaNZ = nHatZ - nHatZM;

   printf("Testing ComputeOutwardUnitNormal: \n");
   printf("nHatX: %12.8f, nHatY: %12.8f, nHatZ: %12.8f\n", nHatX, nHatY, nHatZ);
   printf("Difference between MATLAB and C++:\n");
   printf("deltaNX: %30.20e, deltaNY: %30.20e, deltaNZ: %30.20e\n", deltaNX, deltaNY, deltaNZ);

   return (0);
}
