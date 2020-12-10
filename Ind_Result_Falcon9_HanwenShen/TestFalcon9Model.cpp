#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Falcon9Model.h"
#include "AutopilotModel.h"
#include "AutopilotModel.cpp"
#include "SphericalEarthModel.h"
#include "SphericalEarthModelConstants.h"
#include "Falcon9Model.cpp"


const double DEGTORAD = atan(1.0) / 45.0;
const double RADTODEG = 1.0 / DEGTORAD;

int main(int argc, char **argv)
{
   FILE *output = fopen("Falcon9Model.dat", "w");

   Falcon9Model rocket;

   //   Cape Canaveral Air Force Station Space Launch Complex 40 - Cape Canaveral, FL.

   double launchLatitude  =  28.5620 * DEGTORAD;
   double launchLongitude = -80.5772 * DEGTORAD;
   double launchAltitude  =  32.41;               //   Center of mass of rocket.

   double orbitInclination = 51.6 * DEGTORAD;

   double launchBearing = asin(cos(orbitInclination)/ cos(launchLatitude));

   rocket.Initialize(launchLatitude, launchLongitude, launchAltitude, launchBearing);

   double t      =   0.0;
   double deltaT =   0.05;
   double tEnd   = 160.0;

   //   Output initial state.

   fprintf(output, "%30.20e ", t);

   for (int k = 0; k < 16; k++)
   {
      fprintf(output, "%30.20e ", rocket.u[k]);
   }

   double vXB = 0.0;
   double vYB = 0.0;
   double vZB = 0.0;

   //   Output additional data.

   fprintf(output, "%30.20e %30.20e %30.20e", vXB, vYB, vZB);
   fprintf(output, "%30.20e %30.20e", rocket.deltaPC, rocket.deltaYC);
   fprintf(output, "%30.20e %30.20e %30.20e", rocket.totalAcceleration[0], rocket.totalAcceleration[1], rocket.totalAcceleration[2]);

   fprintf(output, "\n");

   while (t <= tEnd)
   {
      //   Set control parameters at time t.

      //   Set thrust level.

      rocket.eta = rocket.autopilot.ComputeThrustThrottleLevel(t);

      //   Compute velocity in body frame.

      double vX = rocket.u[3];
      double vY = rocket.u[4];
      double vZ = rocket.u[5];

      double qS = rocket.u[12];
      double qX = rocket.u[13];
      double qY = rocket.u[14];
      double qZ = rocket.u[15];

      double xHatB[3] = {0.0};
      double yHatB[3] = {0.0};
      double zHatB[3] = {0.0};

      rocket.ConvertQuaterionToBodyAxes(qS, qX, qY, qZ, xHatB, yHatB, zHatB);

      //   Determine velocity in body frame.

      vXB = vX * xHatB[0] + vY * xHatB[1] + vZ * xHatB[2];
      vYB = vX * yHatB[0] + vY * yHatB[1] + vZ * yHatB[2];
      vZB = vX * zHatB[0] + vY * zHatB[1] + vZ * zHatB[2];

      //   Set thrust gimbal angles.

      //   Update PID controllers.

      double q  = rocket.u[10];                                    //   rad/s.
      double qC = rocket.autopilot.ComputeCommandedPitchRate(t);   //   rad/s.

      if (t < 30.0)
      {
         rocket.deltaPC = rocket.autopilot.ComputeCommandedGimbalPitchRate(q, qC, deltaT);
      }
      else if (t < 160.0)
      {
         rocket.deltaPC = rocket.autopilot.ComputeCommandedGimbalGravityTurnPitch(vZB, deltaT);
         rocket.deltaYC = 0.0;//rocket.autopilot.ComputeCommandedGimbalGravityTurnYaw(vYB, deltaT);
      }
      else
      {
         rocket.deltaPC = 0.0;
         rocket.deltaYC = 0.0;
      }

      //   Update rocket state.

      double tNew = t + deltaT;

      rocket.UpdateState(t, tNew, rocket.u, rocket.u);

      t = tNew;

      //   Output updated state.

      fprintf(output, "%30.20e ", t);
      //printf("testa");
      for (int k = 0; k < 16; k++)
      {
         fprintf(output, "%30.20e ", rocket.u[k]);
      }

      //   Output additional data.

      fprintf(output, "%30.20e %30.20e %30.20e", vXB, vYB, vZB);
      fprintf(output, "%30.20e %30.20e", rocket.deltaPC, rocket.deltaYC);
      fprintf(output, "%30.20e %30.20e %30.20e", rocket.totalAcceleration[0], rocket.totalAcceleration[1], rocket.totalAcceleration[2]);

      fprintf(output, "\n");
   }

   fclose(output);

   return (0);
}
