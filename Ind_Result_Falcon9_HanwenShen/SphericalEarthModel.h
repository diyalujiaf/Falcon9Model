//   SphericalEarthModel.h

#ifndef SPHERICAL_EARTH_MODEL
#define SPHERICAL_EARTH_MODEL

class SphericalEarthModel
{
public:

   SphericalEarthModel();    // Constructor

   ~SphericalEarthModel();   // Destructor

   void ComputeECEFFromLLA(double  latitude,
                           double  longitude,
                           double  altitude,
                           double &xECEF,
                           double &yECEF,
                           double &zECEF);


   void ComputeLLAFromECEF(double  xECEF,
                           double  yECEF,
                           double  zECEF,
                           double &latitude,
                           double &longitude,
                           double &altitude);


   void ConvertENUToECEF(double  enuLatitude,
                         double  enuLongitude,
                         double  enuAltitude,
                         double  xENU,
                         double  yENU,
                         double  zENU,
                         bool    positionData,
                         double &xECEF,
                         double &yECEF,
                         double &zECEF);

   void ComputeGravitationalAcceleration(double  x,
                                         double  y,
                                         double  z,
                                         double &gx,
                                         double &gy,
                                         double &gz);

   void ComputeAtmosphere(double  altitude,            //    Altitude of missile (m)
                          double &airDensity,          //    Air density (kg/(m^3))
                          double &speedOfSound,        //    Speed of sound (m/s)
                          double &pressure,            //    Pressure (N/(m^2))
                          double &temperature,         //    Temperature (K)
                          double &dynamicViscosity);   //    Dynamic viscosity (kg/(ms)).

   void ComputeOutwardUnitNormal(double  latitude,
                                 double  longitude,
                                 double &nHatX,
                                 double &nHatY,
                                 double &nHatZ);

   double GetRE();

   double GetG0();

   double GetSeaLevelPressure();

private:
   //   Add remaining constants at the top of the atmosphere model here.  Change R0 in that function to RE.
   static const double RE;
   static const double G0;
   static const double PRESSURESL;
   static const double R_STAR;
   static const double M0;
   static const double R;
   static const double ALPHA;
   static const double T0;
   static const double SPEED_OF_SOUND_T0;
   static const double BETA;
   static const double S;
};

#endif
