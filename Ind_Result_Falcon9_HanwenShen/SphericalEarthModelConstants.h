//   SphericalEarthModelConstants.h

#ifndef SPHERICAL_EARTH_MODEL_CONSTANTS
#define SPHERICAL_EARTH_MODEL_CONSTANTS

#include "SphericalEarthModel.h"

const double SphericalEarthModel::PRESSURESL = 101325.0;    //   Sea level pressure - N/m^2 (or Pascals).
const double SphericalEarthModel::RE         = 6356766.0;   //   Effective Earth radius - m
const double SphericalEarthModel::G0         = 9.80665;

//   Add remaining constants at the top of the atmosphere model here.  Change R0 in that function to RE.
const double SphericalEarthModel::R_STAR            = 8314.32;          //   Universal gas constant in N m / (kmol K).
const double SphericalEarthModel::M0                =   28.9644;        //   Mean molecular weight at sea level in kg/kmol.
const double SphericalEarthModel::R                 = R_STAR / M0;      //   Specific gas constant in N m /(kg K).
const double SphericalEarthModel::ALPHA             = G0 * M0 / R_STAR;
const double SphericalEarthModel::T0                = 288.15;           //   Standard sea level temperature in K (15 C).
const double SphericalEarthModel::SPEED_OF_SOUND_T0 = 340.3;            //   m/s at T0.
const double SphericalEarthModel::BETA              = 1.458e-6;         //   kg/(msK^(1/2)).
const double SphericalEarthModel::S                 = 110.4;            //   Sutherland's constant (K).

#endif
