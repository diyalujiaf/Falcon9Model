//   Falcon9Model.h

#ifndef FALCON9_MODEL
#define FALCON9_MODEL

#include "SphericalEarthModel.h"
#include "AutopilotModel.h"

class Falcon9Model
{
public:

   AutopilotModel autopilot;

   SphericalEarthModel earth;

   double u[16];
   double eta;
   double deltaPC;
   double deltaYC;

   double totalAcceleration[3];

   Falcon9Model();    // Constructor

   ~Falcon9Model();   // Destructor

   void Initialize(double launchLatitude,
                   double launchLongitude,
                   double launchAltitude,
                   double launchBearing);

   void UpdateState(double  t,
                    double  tNew,
                    double *u,
                    double *uNew);

   void ComputeModelEOM(double  t,
                        double *u,
                        double *dudt);

   double ComputeSeaLevelThrust(double t);

   void ConvertBodyAxesToQuaternion(double *xHatB,   //   Body x-axis (centerline).
                                    double *yHatB,   //   Body y-axis.
                                    double *zHatB,   //   Body z-axiz.
                                    double *q);

   void ConvertQuaterionToBodyAxes(double  qW,
                                   double  qX,
                                   double  qY,
                                   double  qZ,
                                   double *xHatB,
                                   double *yHatB,
                                   double *zHatB);

private:

   double dragonMass;
   double falcon9Mass;
   double initialMass;

   double lHatX;
   double lHatY;
   double lHatZ;

   double G0;
   double ISP;
   double TAU;
   double LTV;
   double LENGTH_OF_ROCKET;
   double IXX;
   double I;
   double NU;
   double PRESSURE_SL;
   double RE;

   double lN;
   double dNC;
   double A;
   double AE;
   double AW;
   double fR;
   double pi;
   double ANT;
   double C;
   double machA;
   double machB;
   double criticalReynoldsNumber;
   double skinFrictionFactor;
   double AB;
   double AP;
   double waveDragReductionFactor;
   double baseToReferenceAreaRatio;
   double planformToreferenceAreaRatio;
   double tData[12];    //   Seconds.

   double thrustData[12];   //   Newtons.

   void ComputeAerodynamicCoefficients(double  mach,        //   Mach number.
                                       double  Re,          //   Reynolds number.
                                       double  alpha,       //   Angle of attack (radians).
                                       bool    powered,     //   1 if rocket engine on, 0 otherwise.
                                       double &cA,          //   Axial force coefficient.
                                       double &cN);         //   Normal force coefficient.

   double ComputeZeroLiftDragCoefficient(double  mach,
                                         double  reynoldsNumber,
                                         bool    powered);

   double ComputeCDWB0(double  mach,
                       double  PO,
                       double  C);

   double ComputeCDC(double machC);

};

#endif
