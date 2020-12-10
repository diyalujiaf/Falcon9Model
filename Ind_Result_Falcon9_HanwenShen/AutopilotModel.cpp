#include <math.h>
#include <stdio.h>

#include "AutopilotModel.h"
#include "SphericalEarthModel.h"
#include "SphericalEarthModelConstants.h"


AutopilotModel::AutopilotModel()   // Constructor
{
   oldErrorP   = 0.0;
   errorP      = 0.0;
   integralP   = 0.0;
   derivativeP = 0.0;

   oldErrorY   = 0.0;
   errorY      = 0.0;
   integralY   = 0.0;
   derivativeY = 0.0;

   //   Pitch rate PID Gains.

   KPPR = 5.95;          //   Seconds.
   TIPR = 0.2;
   TDPR = 0.1;
   KIPR = KPPR / TIPR;   //   Dimensionless.
   KDPR = KPPR * TDPR;   //   Seconds^2.

   //   Gravity turn PID Gains.

   KPGT = 0.002;         //   rad s / m.
   TIGT = 2.0;
   TDGT = 1.0;
   KIGT = KPGT / TIGT;   //   rad / m.
   KDGT = KPGT * TDGT;   //   rad s^2 /m.
};

AutopilotModel::~AutopilotModel()   // Destructor
{
};

double AutopilotModel::ComputeCommandedGimbalPitchRate(double  q,        //   Achieved pitch rate (rad/s).
                                                       double  qC,       //   Commanded pitch rate (rad/s).
                                                       double  deltaT)   //   t step.
{
   errorP = qC - q;   //   Commanded - Achieved (rad/s).

   double dError = errorP - oldErrorP;

   oldErrorP = errorP;

   derivativeP = dError / deltaT;
   integralP   = integralP + errorP * deltaT;

   double deltaPC = KPPR * errorP + KIPR * integralP + KDPR * derivativeP;   //   Radians.

   return deltaPC;
};

double AutopilotModel::ComputeCommandedGimbalGravityTurnPitch(double  vZB,      //   Velocity along z-axis.
                                                              double  deltaT)   //   t step.
{
   errorP = -vZB;   //   Commanded (zero) - Achieved (m/s).

   double dError = errorP - oldErrorP;

   oldErrorP = errorP;

   derivativeP = dError / deltaT;
   integralP   = integralP + errorP * deltaT;

   double deltaPC = KPGT * errorP + KIGT * integralP + KDGT * derivativeP;   //   Radians.

   return deltaPC;
};

double AutopilotModel::ComputeCommandedGimbalGravityTurnYaw(double  vYB,      //   Velocity along y-axis.
                                                            double  deltaT)   //   t step.
{
   errorY = -vYB;   //   Commanded (zero) - Achieved (m/s).

   double dError = errorY - oldErrorY;

   oldErrorY = errorY;

   derivativeY = dError / deltaT;
   integralY   = integralY + errorY * deltaT;

   double deltaYC = KPGT * errorY + KIGT * integralY + KDGT * derivativeY;   //   Radians.

   return deltaYC;
};

double AutopilotModel::ComputeThrustThrottleLevel(double t)
{
   double tT[5]   = {0.0, 40.0, 45.0,  68.0, 73.0};
   double etaT[5] = {1.0,  1.0,  0.75,  0.75, 1.0};

   double eta = 1.0;

   if (t < tT[4])
   {
      int k = 0;

      while ((t >= tT[k]) && (k < 4))
      {
         k++;
      }

      double deltaT   = tT[k] - tT[k - 1];
      double deltaEta = etaT[k] - etaT[k - 1];

      eta = etaT[k - 1] + (t - tT[k - 1]) * deltaEta / deltaT;
   }
   return eta;
}

double AutopilotModel::ComputeCommandedPitchRate(double t)
{
   double qC = 0.0;

   if ((t >= 14.0) && (t <= 30.0))
   {
      qC = -0.005454154;
   }

   double factor = 1.5;

   qC = qC * factor;

   return qC;
};

