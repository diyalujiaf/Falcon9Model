#include <math.h>
#include <stdio.h>

#include "SphericalEarthModel.h"
#include "SphericalEarthModelConstants.h"


SphericalEarthModel::SphericalEarthModel()   // Constructor
{
};

SphericalEarthModel::~SphericalEarthModel()   // Destructor
{
};

void SphericalEarthModel::ComputeECEFFromLLA(double  latitude,
                                             double  longitude,
                                             double  altitude,
                                             double &xECEF,
                                             double &yECEF,
                                             double &zECEF)
{
    double range;
    double cosLat;

range = RE + altitude;

cosLat = cos(latitude);

xECEF = range * cosLat * cos(longitude);
yECEF = range * cosLat * sin(longitude);
zECEF = range * sin(latitude);
};

void SphericalEarthModel::ComputeLLAFromECEF(double  x,
                                             double  y,
                                             double  z,
                                             double &latitude,
                                             double &longitude,
                                             double &altitude)
{
    double horizontalRangeSquared;
    double horizontalRange;
    double range;

horizontalRangeSquared = x * x + y * y;
horizontalRange        = sqrt(horizontalRangeSquared);
range                  = sqrt(horizontalRangeSquared + z * z);
altitude               = range - RE;

longitude = atan2(y, x);
latitude  = atan2(z, horizontalRange);
};

void SphericalEarthModel::ConvertENUToECEF(double  enuLatitude,
                                           double  enuLongitude,
                                           double  enuAltitude,
                                           double  xENU,
                                           double  yENU,
                                           double  zENU,
                                           bool    positionData,
                                           double &xECEF,
                                           double &yECEF,
                                           double &zECEF)
{
double cosLat;
double sinLat;
double cosLon;
double sinLon;
cosLat = cos(enuLatitude);
sinLat = sin(enuLatitude);
cosLon = cos(enuLongitude);
sinLon = sin(enuLongitude);

if (positionData)
   {
   zENU = zENU + enuAltitude;

   xECEF = (zENU + RE) * cosLat * cosLon - yENU * sinLat * cosLon - xENU * sinLon;
   yECEF = (zENU + RE) * cosLat * sinLon - yENU * sinLat * sinLon + xENU * cosLon;
   zECEF = (zENU + RE) * sinLat + yENU * cosLat;
   }
else
   {
   xECEF = (zENU * cosLat - yENU * sinLat) * cosLon - xENU * sinLon;
   yECEF = (zENU * cosLat - yENU * sinLat) * sinLon + xENU * cosLon;
   zECEF =  zENU * sinLat + yENU * cosLat;
   }
};

void SphericalEarthModel::ComputeGravitationalAcceleration(double  x,
                                                           double  y,
                                                           double  z,
                                                           double &gx,
                                                           double &gy,
                                                           double &gz)
{
// G0  = 9.80665;
    double earthCenteredRange;
    double altitude;
    double factor;
    double gravityAtAltitude;
    double xHat;

    earthCenteredRange = sqrt(x * x + y * y + z * z);
    altitude           = earthCenteredRange - RE;

    factor            = 1.0 + altitude / RE;
    gravityAtAltitude = G0 / (factor * factor);
    gx                = -gravityAtAltitude * x / earthCenteredRange;
    gy                = -gravityAtAltitude * y / earthCenteredRange;
    gz                = -gravityAtAltitude * z / earthCenteredRange;

};


void SphericalEarthModel::ComputeAtmosphere(double  altitude,           //    Altitude of missile (m)
                                            double &airDensity,         //    Air density (kg/(m^3))
                                            double &speedOfSound,       //    Speed of sound (m/s)
                                            double &pressure,           //    Pressure (N/(m^2))
                                            double &temperature,        //    Temperature (K)
                                            double &dynamicViscosity)   //    Dynamic viscosity (kg/(ms))
{

    double geoAltitude = altitude / (1.0 + altitude / RE);
    double lnPressure;

if (geoAltitude <= 11000.0) //   Troposphere.
   {
    temperature = T0 - 0.0065 * geoAltitude;
    lnPressure  = -ALPHA * log(T0 / temperature) / 0.0065;
    pressure    = PRESSURESL * exp(lnPressure);
   }

else if (geoAltitude <= 20000.0) //  Tropsphere.
   {
   temperature = 216.65;
   pressure    = 22664.332136278 * exp(-ALPHA * (geoAltitude - 11000.0) / temperature);
    }
else if (geoAltitude <= 32000.0) //   Stratosphere.
   {
   temperature = 216.65 + 0.001 * (geoAltitude - 20000.0);
   lnPressure  = ALPHA * log(216.65 / temperature) / 0.001;
   pressure    = 5485.3789061267 * exp(lnPressure);
   }
else if (geoAltitude <= 47000.0) //   Stratosphere.
   {
   temperature = 228.65 + 0.0028 * (geoAltitude - 32000.0);
   lnPressure  = ALPHA * log(228.65 / temperature) / 0.0028;
   pressure    = 869.92672346467 * exp(lnPressure);
   }
else if (geoAltitude <= 51000.0) //   Stratosphere.
   {
   temperature = 270.65;
   pressure    = 111.15136951935 * exp(-ALPHA * (geoAltitude - 47000.0) / temperature);
   }
else if (geoAltitude <= 71000.0) //  Mesosphere.
   {
   temperature = 270.65 - 0.0028 * (geoAltitude - 51000.0);
   lnPressure  = -ALPHA * log(270.65 / temperature) / 0.0028;
   pressure    = 67.107449909353 * exp(lnPressure);
   }
else if (geoAltitude <= 84852.0) //  Mesosphere.
   {
   temperature = 214.65 - 0.002 * (geoAltitude - 71000.0);
   lnPressure  = -ALPHA * log(214.65 / temperature) / 0.002;
   pressure    = 3.9676007402804 * exp(lnPressure);
   }
else
   {
   temperature = 186.95;
   pressure    = 0.37510225611954 * exp(-ALPHA * (geoAltitude - 84852.0) / temperature);
   }
   airDensity       = pressure / (R * temperature);
speedOfSound     = SPEED_OF_SOUND_T0 * sqrt(temperature / T0);
dynamicViscosity = BETA * pow(temperature,1.5) / (temperature + S);

//(temperature^1.5)
};

void SphericalEarthModel::ComputeOutwardUnitNormal(double  latitude,
                                                   double  longitude,
                                                   double &nHatX,
                                                   double &nHatY,
                                                   double &nHatZ)
{
   //   This function computes a unit vector for the local vertical.
       double cosLat;
cosLat = cos(latitude);

nHatX = cosLat * cos(longitude);
nHatY = cosLat * sin(longitude);
nHatZ = sin(latitude);
};


double SphericalEarthModel::GetRE()
{
   return RE;
};

double SphericalEarthModel::GetG0()
{
   return G0;
};

double SphericalEarthModel::GetSeaLevelPressure()
{
   return PRESSURESL;
};

//double SphericalEarthModel::GetR_STAR()
//{
//   return R_STAR;
//};
//
//double SphericalEarthModel::GetM0()
//{
//   return M0;
//};
//
//double SphericalEarthModel::GetR()
//{
//   return R;
//};
//
//double SphericalEarthModel::GetALPHA()
//{
//   return ALPHA;
//};
//
//double SphericalEarthModel::GetT0()
//{
//   return T0;
//};
//
//double SphericalEarthModel::GetSPEED_OF_SOUND_T0()
//{
//   return SPEED_OF_SOUND_T0;
//};
//
//double SphericalEarthModel::GetBETA()
//{
//   return BETA;
//};
//
//double SphericalEarthModel::GetS()
//{
//   return S;
//};

