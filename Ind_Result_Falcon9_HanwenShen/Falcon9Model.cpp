#include <math.h>
#include <stdio.h>

#include "Falcon9Model.h"
#include "AutopilotModel.h"
#include "SphericalEarthModel.h"
#include "SphericalEarthModelConstants.h"
#include "SphericalEarthModel.cpp"

#include <algorithm>


Falcon9Model::Falcon9Model()   // Constructor
{
   //   Initialize masses.

   dragonMass  =   6000.0;
   falcon9Mass = 549054.0;

   initialMass  = falcon9Mass + dragonMass;

   //   Initialize thrust throttle.

   eta = 1.0;

   //   Initialize thrust gimbal angles.

   deltaPC = 0.0;
   deltaYC = 0.0;

   lHatX = 0.0;
   lHatY = 0.0;
   lHatZ = 0.0;

   totalAcceleration[0] = 0.0;
   totalAcceleration[1] = 0.0;
   totalAcceleration[2] = 0.0;

   AE  = 5.526;   //   Nozzle exit area (m^2).
   ISP = 282.0;   //   Sea level specific impulse.
   TAU = 0.25;    //   TVC time constant.
   LTV = 32.41;   //   Distance from center of mass to end of rocket (m).

   LENGTH_OF_ROCKET = 64.82;   //  Meters.  This is for the full rocket.  It needs to change per stage.

   //   These are just estimated from a cylinder

   IXX =    919363.47;   //   Roll moment of inertia.
   I   = 147140574.24;   //   Pitch/Yaw moment of inertia (I = IYY = IZZ).
   NU  = IXX / I;

   //   Constants needed from spherical Earth model.

   PRESSURE_SL = earth.GetSeaLevelPressure();
   RE          = earth.GetRE();
   G0          = earth.GetG0();

   //   Constants used in aerodynamic coefficient model.

   lN  =   6.88;     //   Length of nose cone (m)
   dNC =   4.0;      //   Diameter of nose cone (m)
   A   =  10.52;     //   Aerodynamic reference area (m^2)
   AW  = 651.03;     //   Total wetted area for entire rocket during ascent phase (m^2).
   fR  = lN / dNC;   //   Fineness ratio (for sharp cone).

   pi = 4.0 * atan(1.0);

   ANT = (pi * dNC * dNC) / 4.0; // Area of the nose tip (m^2)

   C = (((A - ANT) / A) * pow((atan(0.5 / fR)), 1.69) + (ANT / A) * 0.665); // constant

   machA = 0.8;
   machB = 1.2;

   criticalReynoldsNumber = 5.0e5;
   skinFrictionFactor     = (AW / A) * 0.074;

   AB =  10.52;   //   Cross-sectional area of base (m^2).
   AP = 232.37;   //   Planform area for entire rocket during ascent phase (m^2).

   waveDragReductionFactor      = AE / A;
   baseToReferenceAreaRatio     = AB / A;
   planformToreferenceAreaRatio = AP / A;


};

Falcon9Model::~Falcon9Model()   // Destructor
{
};

void Falcon9Model::Initialize(double launchLatitude,
                              double launchLongitude,
                              double launchAltitude,
                              double launchBearing)
{
   double initialX = 0.0;
   double initialY = 0.0;
   double initialZ = 0.0;

   earth.ComputeECEFFromLLA(launchLatitude, launchLongitude, launchAltitude, initialX, initialY, initialZ);

   earth.ComputeOutwardUnitNormal(launchLatitude, launchLongitude, lHatX, lHatY, lHatZ);

   double cosLaunchBearing = cos(launchBearing);
   double sinLaunchBearing = sin(launchBearing);

   double lHatPerpX = 0.0;
   double lHatPerpY = 0.0;
   double lHatPerpZ = 0.0;

   earth.ConvertENUToECEF(launchLatitude, launchLongitude, 0.0, sinLaunchBearing, cosLaunchBearing, 0.0, false, lHatPerpX, lHatPerpY, lHatPerpZ);

   double lHatPerp[3] = {lHatPerpX, lHatPerpY, lHatPerpZ};

   //   Initial orientation.

   double xHatB0[3] = {lHatX, lHatY, lHatZ};
   double yHatB0[3] = {lHatPerpY * lHatZ - lHatPerpZ * lHatY, lHatPerpZ * lHatX - lHatPerpX * lHatZ, lHatPerpX * lHatY - lHatPerpY * lHatX};
   double zHatB0[3] = {lHatPerpX, lHatPerpY, lHatPerpZ};

   printf("%12.8f %12.8f %12.8f\n", xHatB0[0], xHatB0[1], xHatB0[2]);
   printf("%12.8f %12.8f %12.8f\n", yHatB0[0], yHatB0[1], yHatB0[2]);
   printf("%12.8f %12.8f %12.8f\n", zHatB0[0], zHatB0[1], zHatB0[2]);

   double q0[4] = {0.0};

   ConvertBodyAxesToQuaternion(xHatB0, yHatB0, zHatB0, q0);

   //   State vector: pos, vel, mass, deltas, omega, quaternion.

   u[0]  = initialX;      //   Position (m).
   u[1]  = initialY;
   u[2]  = initialZ;
   u[3]  = 0.0;           //   Velocity (m/s).
   u[4]  = 0.0;
   u[5]  = 0.0;
   u[6]  = initialMass;   //   Mass (kg).
   u[7]  = 0.0;           //   Thrust gimbal angles (rad).
   u[8]  = 0.0;
   u[9]  = 0.0;           //   Angular velocity (rad/s).
   u[10] = 0.0;
   u[11] = 0.0;
   u[12] = q0[0];         //   Orientation as unit quaternion.
   u[13] = q0[1];
   u[14] = q0[2];
   u[15] = q0[3];

};

void Falcon9Model::UpdateState(double  t,
                               double  tNew,
                               double *u,
                               double *uNew)
{
   //   This is a 4th-order Runge-Kutta method.  It is set up for this specific case of a state
   //   vector of length 16.

   int n = 16;

   double f[16];
   double f1[16];
   double f2[16];
   double f3[16];
   double f4[16];
   double u1[16];

   double h = tNew - t;

   ComputeModelEOM(t, u, f);

   for (int k = 0; k < n; k++)
   {
       //printf("%f",12);
      f1[k] = h * f[k];

      u1[k] = u[k] + 0.5 * f1[k];
   }

   t = t + (double)0.5 * h;

   ComputeModelEOM(t, u1, f);

   for (int k = 0; k < n; k++)
   {
      f2[k] = h * f[k];
      u1[k] = u[k] + (double)0.5 * f2[k];
   }

   ComputeModelEOM(t, u1, f);

   for (int k = 0; k < n; k++)
   {
      f3[k] = h * f[k];
      u1[k] = u[k] + f3[k];
   }

   t = tNew;

   ComputeModelEOM(t, u1, f);

   for (int k = 0; k < n; k++)
   {
      f4[k] = h * f[k];

      //  Update state.

      uNew[k] = u[k] + (f1[k] + 2.0 * (f2[k] + f3[k]) + f4[k]) / 6.0;
   }

   //   Save total acceleration from last EOM call.

   totalAcceleration[0] = f[3];
   totalAcceleration[1] = f[4];
   totalAcceleration[2] = f[5];
};


void Falcon9Model::ComputeModelEOM(double  t,
                                   double *u,
                                   double *dudt)
{



    double AE  = 5.526;   //   Nozzle exit area (m^2).
    double G0  = 9.80665;
    double ISP = 282.0;   //   Sea level specific impulse.
    double TAU = 0.25;    //   TVC time constant.
    double LTV = 32.41;   //   Distance from center of mass to end of rocket (m).

//   These are just estimated from a cylinder

    double IXX	=    919363.47;   //
    double I   = 147140574.24;   //  I = IYY = IZZ;
    double NU  = IXX / I;

    double referenceArea = 10.52;   //   m^2.

    double PRESSURE_SL = 101325.0;

   double x      = u[0];
   double y      = u[1];
   double z      = u[2];
   double vX     = u[3];
   double vY     = u[4];
   double vZ     = u[5];
   double mass   = u[6];
   double deltaP = u[7];
   double deltaY = u[8];
   double omegaX = u[9];
   double omegaY = u[10];
   double omegaZ = u[11];
   double qS     = u[12];
   double qX     = u[13];
   double qY     = u[14];
   double qZ     = u[15];


   //   ADD THIS REST OF THE MODEL CODE HERE.
    double mDot;
    double thrustSL;
    double pressureThrust;
    double airDensity;
    double pressure;
    double temperature;
    double thrustAtAltitude;
    double speed;
    double mach;
    double reynoldsNumber;
    double alpha;
    bool powered;
    double cA;
    double cN;
    double axialForce;
    double normalForce;
    double dynamicPressure;
    double speedOfSound;
    double dynamicViscosity;
    double cosAlpha;
    double axialThrust;
    double normalAccMag;
    double axialAcc[3];
    double verticalAcceleration;
    double normalAcc[3];
    double axialAccMag;
    double gX;
    double gY;
    double gZ;


    double cosDeltaP = cos(deltaP);
    double sinDeltaP = sin(deltaP);
    double cosDeltaY = cos(deltaY);
    double sinDeltaY = sin(deltaY);
    double xHatB[3];
    double yHatB[3];
    double zHatB[3];

    ConvertQuaterionToBodyAxes(qS, qX, qY, qZ,xHatB,yHatB,zHatB);


    double earthCenteredRange = sqrt(x * x + y * y + z * z);
    double speedSquared       = vX * vX + vY * vY + vZ * vZ;
    double altitude           = earthCenteredRange - RE;


    double totalAccelerationX= totalAcceleration[0];
    double totalAccelerationY= totalAcceleration[1];
    double totalAccelerationZ= totalAcceleration[2];
    //   Compute sea level thrust.
    earth.ComputeAtmosphere(altitude,airDensity, speedOfSound, pressure, temperature, dynamicViscosity);

//printf("test%f,",pressure);
    thrustSL = ComputeSeaLevelThrust(t);

//   Throttle thrust by eta.

    thrustSL = eta * thrustSL;

//   Compute mass flow rate.

    mDot = thrustSL / (G0 * ISP);

    thrustAtAltitude = 0.0;


    if (thrustSL > 0.0)
    {
        pressureThrust   = (PRESSURE_SL - pressure) * AE;
        thrustAtAltitude = thrustSL + pressureThrust;

    }
//   Compute aerodynamic forces.
    dynamicPressure = 0.5 * airDensity * speedSquared;

    speed = sqrt(speedSquared);
    mach  = speed / speedOfSound;

//   Compute Reynolds number.

    reynoldsNumber = airDensity * speed * LENGTH_OF_ROCKET / dynamicViscosity;

//   Compute total angle of attack.

    alpha = 0.0;

    if (speed > 0.0)
    {
        cosAlpha = (vX * xHatB[0] + vY * xHatB[1] + vZ * xHatB[2]) / speed;

   //   Fix round-off error.

        cosAlpha = std::min(std::max(cosAlpha, -1.0), 1.0);

        alpha = acos(cosAlpha);
    }

//powered = thrustAtAltitude > 0.0;

    if(thrustAtAltitude > 0.0)
    {
        powered = thrustAtAltitude ;
    }
    else
    {
        powered = false;
    }


    ComputeAerodynamicCoefficients(mach, reynoldsNumber, alpha, powered,cA, cN);


    axialForce  = dynamicPressure * cA * referenceArea;
    normalForce = dynamicPressure * cN * referenceArea;


    earth.ComputeGravitationalAcceleration(x, y, z,gX,gY,gZ);

//   Axial thrust.

    axialThrust = thrustAtAltitude * cosDeltaP * cosDeltaY;



    axialAccMag  =  (axialThrust - axialForce) / mass;
    normalAccMag =  normalForce / mass;
//from array to list
    axialAcc[0] =  axialAccMag * xHatB[0];
    axialAcc[1] =  axialAccMag * xHatB[1];   //   This should work fine for this application.
    axialAcc[2] =  axialAccMag * xHatB[2];
//   This needs to change sign with pitch angle of attack.

    normalAcc[0] = - normalAccMag * zHatB[0];
    normalAcc[1] = - normalAccMag * zHatB[1];   //   This should work fine for this application.
    normalAcc[2] = - normalAccMag * zHatB[2];

    double v[3];
    v[0]=vX;
    v[1]=vY;
    v[2]=vZ;

    if ((v[0]*zHatB[0]+ v[1]*zHatB[1] + v[2]*zHatB[2]) < 0.0)   //   Negative pitch angle of attack
    {
   //   Direction of normal acceleration in positive z axis.

        normalAcc[0] = - normalAcc[0];
        normalAcc[1] = - normalAcc[1];   //   This should work fine for this application.
        normalAcc[2] = - normalAcc[2];
    }
    double gravityAcc[3];
    gravityAcc[0] = gX;
    gravityAcc[1] = gY;
    gravityAcc[2] = gZ;
//fprintf('%12.3f %12.3f %12.3f\n', t, axialAccMag, thrustAtAltitude);

    totalAccelerationX = axialAcc[0] + normalAcc[0] + gravityAcc[0];
    totalAccelerationY = axialAcc[1] + normalAcc[1] + gravityAcc[1];
    totalAccelerationZ = axialAcc[2] + normalAcc[2] + gravityAcc[2];
    totalAcceleration[0] = totalAccelerationX;
    totalAcceleration[1] = totalAccelerationY;
    totalAcceleration[2] = totalAccelerationZ;

    verticalAcceleration = totalAccelerationX*lHatX + totalAccelerationY*lHatY  + totalAccelerationZ*lHatZ ;

//fprintf('%8.3f %8.3f\n', t, verticalAcceleration);

    if ((t < 5.0) && (verticalAcceleration < 0.0))
    {
        totalAccelerationX= 0.0;
        totalAccelerationY= 0.0;
        totalAccelerationZ= 0.0;
        totalAcceleration[0] = totalAccelerationX;
        totalAcceleration[1] = totalAccelerationY;
        totalAcceleration[2] = totalAccelerationZ;
    }
//   Compute torques.
   //   Compute quaternion cost.

    double tauX = 0.0;//me
    double tauY = LTV * thrustAtAltitude * sinDeltaP;//me
    double tauZ = LTV * thrustAtAltitude * cosDeltaP * sinDeltaY;//me

   double lambda = 1.0 - (qS * qS + qX * qX + qY * qY + qZ * qZ);

   dudt[0]  =  vX;
   dudt[1]  =  vY;
   dudt[2]  =  vZ;
   dudt[3]  =  totalAccelerationX;
   dudt[4]  =  totalAccelerationY;
   dudt[5]  =  totalAccelerationZ;
   dudt[6]  = -mDot;
   dudt[7]  = (deltaPC - deltaP) / TAU;
   dudt[8]  = (deltaYC - deltaY) / TAU;
   dudt[9]  = tauX / (NU * I);
   dudt[10] = tauY / I + (1.0 - NU) * omegaX * omegaZ;
   dudt[11] = tauZ / I + (1.0 - NU) * omegaX * omegaY;
   dudt[12] = -0.5 * (omegaX * qX + omegaY * qY + omegaZ * qZ) + lambda * qS;
   dudt[13] =  0.5 * (omegaX * qS + omegaZ * qY - omegaY * qZ) + lambda * qX;
   dudt[14] =  0.5 * (omegaY * qS - omegaZ * qX + omegaX * qZ) + lambda * qY;
   dudt[15] =  0.5 * (omegaZ * qS + omegaY * qX - omegaX * qY) + lambda * qZ;
};

double Falcon9Model::ComputeSeaLevelThrust(double t)
{
    double tData[12] = {0.0, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 152.0, 154.0};    //   Seconds.

    double thrustData[12] = {      0.0, 5437000.0, 5728000.0, 6030000.0, 6296000.0, 6518000.0, 6911000.0,
                            7269000.0, 7529000.0, 7607000.0, 7607000.0,       0.0};   //   Newtons.

   double seaLevelThrust = 0.0;

   if (t < tData[11])
   {
      int k = 0;

      while ((t >= tData[k]) && (k < 11))
      {
         k++;
      }

      double deltaT      = tData[k] - tData[k - 1];
      double deltaThrust = thrustData[k] - thrustData[k - 1];

      seaLevelThrust = thrustData[k - 1] + (t - tData[k - 1]) * deltaThrust / deltaT;
   }

   return seaLevelThrust;
};

void Falcon9Model::ConvertBodyAxesToQuaternion(double *xHatB,   //   Body x-axis (centerline).
                                               double *yHatB,   //   Body y-axis.
                                               double *zHatB,   //   Body z-axiz.
                                               double *q)
{
   //   This function converts the body axes to a unit quaterion.

   double EPSILON = 0.000000000001;

   double C[3][3] = {0.0};

   C[0][0] = xHatB[0];
   C[1][0] = xHatB[1];
   C[2][0] = xHatB[2];

   C[0][1] = yHatB[0];
   C[1][1] = yHatB[1];
   C[2][1] = yHatB[2];

   C[0][2] = zHatB[0];
   C[1][2] = zHatB[1];
   C[2][2] = zHatB[2];

   double sSquared = 0.25 * (1.0 + C[0][0] + C[1][1] + C[2][2]);

   if (sSquared > EPSILON)
   {
      q[0] = sqrt(sSquared);
      q[1] = (C[2][1] - C[1][2]) / (4.0 * q[0]);
      q[2] = (C[0][2] - C[2][0]) / (4.0 * q[0]);
      q[3] = (C[1][0] - C[0][1]) / (4.0 * q[0]);
   }
   else
   {
      q[0] = 0.0;

      double xSquared = -0.5 * (C[1][1] + C[2][2]);

      if (xSquared > EPSILON)
	  {
         q[1] = sqrt(xSquared);
         q[2] = C[1][0] / (2.0 * q[1]);
         q[3] = C[2][0] / (2.0 * q[1]);
	  }
      else
	  {
         q[1] = 0.0;

         double ySquared = 0.5 * (1.0 - C[2][2]);

         if (ySquared > EPSILON)
		 {
            q[2] = sqrt(ySquared);
            q[3] = C[2][1] / (2.0 * q[2]);
		 }
         else
		 {
            q[2] = 0.0;
            q[3] = 1.0;
		 }
      }
   }
};

void Falcon9Model::ConvertQuaterionToBodyAxes(double  qW,
                                              double  qX,
                                              double  qY,
                                              double  qZ,
                                              double *xHatB,
                                              double *yHatB,
                                              double *zHatB)
{
   //   This function converts a unit quaternion that represents orientation to
   //   the body axes of the rocket.

   //   Body axes are the columns of the direction cosine matrix.

   xHatB[0] = 1.0 - 2.0 * (qY * qY + qZ * qZ);
   xHatB[1] =       2.0 * (qX * qY + qW * qZ);
   xHatB[2] =       2.0 * (qX * qZ - qW * qY);

   yHatB[0] =       2.0 * (qX * qY - qW * qZ);
   yHatB[1] = 1.0 - 2.0 * (qX * qX + qZ * qZ);
   yHatB[2] =       2.0 * (qY * qZ + qW * qX);

   zHatB[0] =       2.0 * (qX * qZ + qW * qY);
   zHatB[1] =       2.0 * (qY * qZ - qW * qX);
   zHatB[2] = 1.0 - 2.0 * (qX * qX + qY * qY);
};

void Falcon9Model::ComputeAerodynamicCoefficients(double  mach,        //   Mach number.
                                                  double  Re,          //   Reynolds number.
                                                  double  alpha,       //   Angle of attack (radians).
                                                  bool    powered,     //   1 if rocket engine on, 0 otherwise.
                                                  double &cA,          //   Axial force coefficient.
                                                  double &cN)          //   Normal force coefficient.
{
   double cD0 = ComputeZeroLiftDragCoefficient(mach, Re, powered);

   //   Compute cross flow Mach number.

   double machC = mach * sin(alpha);

   //   Compute crossflow drag coefficient for circular cylinder.

   double cDc = ComputeCDC(machC);

   //   Modified cD0 for angle of attack to get cA.

   //   Drag proportionality factor (based on length-to-diameter ratio of 16.6).

   double eta0 = 0.739;   //   See Allen 1949.

   double eta = 1.0;

   if (machC <= 1.8)
   {
      eta = ((1.0 - eta0) / 1.8) * machC + eta0;
   }

   double factorA = baseToReferenceAreaRatio * alpha;
   double factorB = eta * cDc * planformToreferenceAreaRatio * alpha * alpha;

   cA = cD0 + (factorA + factorB) * alpha;

   //   Compute normal force coefficient.

   cN = 2.0 * factorA + factorB;
};

double Falcon9Model::ComputeZeroLiftDragCoefficient(double  mach,
                                                    double  reynoldsNumber,
                                                    bool    powered)
{
   double cD0   = 0.0;
   double cDWB0 = 0.0;
   double PO    = 1.0;

   if (powered)
   {
      PO = PO - waveDragReductionFactor;
   }

   if ((mach <= machA) || (mach >= machB))
   {
      //   Equations for wave and base drag coefficients.

      cDWB0 = ComputeCDWB0(mach, PO, C);
   }
   else   //   Transonic region.
   {
      //   Evaluate wave and base drag coeffcients at machA and machB.

      double cDWB0A = ComputeCDWB0(machA, PO, C);
      double cDWB0B = ComputeCDWB0(machB, PO, C);

      //   Evaluate wave and base drag coeffcients at machA and machB.

      double dCDWB0dMA = PO * 0.26 * machA;
      double dCDWB0dMB = -(C * 3.668 / (machB * machB * machB) + PO * 0.25 / (machB * machB));

      //   Cubic Hermite interpolation.

      double deltaMachA         = mach  - machA;
      double deltaMachASquared  = deltaMachA * deltaMachA;
      double deltaMachB         = machB - mach;
      double deltaMachBSquared  = deltaMachB * deltaMachB;
      double deltaMachAB        = machB - machA;
      double deltaMachABSquared = deltaMachAB * deltaMachAB;

      double factor0 = (1.0 + 2.0 * deltaMachA / deltaMachAB) * deltaMachBSquared / deltaMachABSquared;
      double factor1 = (1.0 + 2.0 * deltaMachB / deltaMachAB) * deltaMachASquared / deltaMachABSquared;

      double factorP0 =  deltaMachA * deltaMachBSquared / deltaMachABSquared;
      double factorP1 = -deltaMachASquared * deltaMachB / deltaMachABSquared;

      cDWB0 = factor0  * cDWB0A + factor1  * cDWB0B + factorP0 * dCDWB0dMA + factorP1 * dCDWB0dMB;
   }

   //   Compute skin friction coefficient.

   if (reynoldsNumber < criticalReynoldsNumber)
   {
      reynoldsNumber = criticalReynoldsNumber;
   }

   double cDF = skinFrictionFactor / pow(reynoldsNumber, 0.2);

   cD0 = cDWB0 + cDF;

   return cD0;
};

double Falcon9Model::ComputeCDWB0(double  mach,
                                  double  PO,
                                  double  C)
{
   //   This function compute the combined wave and base drag coefficients
   //  for the zero-lift drag coefficient.

   double cDWB0 = 0.0;

   if (mach < 1.0)
   {
      cDWB0 = PO * (0.12 + 0.13 * mach * mach);
   }
   else
   {
      cDWB0 = C * (1.586 + 1.834 / (mach * mach)) + PO * 0.25 / mach;
   }

   return cDWB0;
};

double Falcon9Model::ComputeCDC(double machC)   //   Crossflow Mach number.
{
   //   This computes the crossflow drag coefficient for a circular cylinder
   //   based on Allen 1949 Figure 9.

   double machCT[18] = {0.0,    0.3197, 0.3813, 0.4939, 0.6085, 0.6777, 0.7325, 0.7755, 0.8256,
                        0.8733, 0.9783, 1.148,  1.267,  1.489,  1.67,   1.854,  1.988,  2.202};

   double cDcT[18]  = {1.2,   1.2,   1.2249, 1.3252, 1.526, 1.686, 1.768, 1.8,   1.818,
                       1.822, 1.795, 1.698,  1.61,   1.463, 1.371, 1.31,  1.273, 1.227};

   double cDc = 1.2;

   if (machC < machCT[17])
   {
      int i = 0;

      while ((machC >= machCT[i]) && (i < 17))
      {
         i++;
      }

      double deltaMachC = machCT[i] - machCT[i - 1];
      double deltaCDC   = cDcT[i] - cDcT[i - 1];

      cDc = cDcT[i - 1] + (machC - machCT[i - 1]) * deltaCDC / deltaMachC;
   }

   return cDc;
};
