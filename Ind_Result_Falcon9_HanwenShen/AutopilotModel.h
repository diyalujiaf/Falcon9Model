//   AutopilotModel.h

#ifndef AUTOPILOT_MODEL
#define AUTOPILOT_MODEL

class AutopilotModel
{
public: 

   AutopilotModel();    // Constructor

   ~AutopilotModel();   // Destructor
   
   double ComputeCommandedGimbalPitchRate(double  q,         //   Achieved pitch rate (rad/s).
                                          double  qC,        //   Commanded pitch rate (rad/s).
                                          double  deltaT);   //   t step.


   double ComputeCommandedGimbalGravityTurnPitch(double  vZB,     //   Velocity along z-axis.
                                               double  deltaT);   //   t step.

   double ComputeCommandedGimbalGravityTurnYaw(double  vYB,       //   Velocity along y-axis.
                                               double  deltaT);   //   t step.

   double ComputeThrustThrottleLevel(double t);

   double ComputeCommandedPitchRate(double t);


private:

   double oldErrorP;
   double errorP;
   double integralP;
   double derivativeP;
   
   double oldErrorY;
   double errorY;
   double integralY;
   double derivativeY;

   //   Pitch rate PID gain parameters.

   double KPPR;
   double TIPR;
   double TDPR;
   double KIPR;
   double KDPR;

   //   Gravity turn PID gain parameters.

   double KPGT;
   double TIGT;
   double TDGT;
   double KIGT;
   double KDGT;
};

#endif
