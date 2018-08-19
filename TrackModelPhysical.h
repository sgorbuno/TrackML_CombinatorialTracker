#ifndef TRACKMODELPHYSICAL_H
#define TRACKMODELPHYSICAL_H

#include "Geo.h"

struct HitMC;

class TrackModelPhysical
{
 public:
  double x0; // position of the parameterization point
  double y0;
  double z0;
  double ex; // directions of u axis
  double ey;
  double pu; // momentum
  double pv;
  double pz; 
  double  q;
  double pt;  
  double ptInv;
   
  //

  int createXY(double xf, double yf, double zf,
	       double xm, double ym, double zm,
	       double xl, double yl, double zl, double Bz  );

  int createXYFirst(double xf, double yf, double zf,
		    double xm, double ym, double zm,
		    double xl, double yl, double zl, double Bz  );

  int createXYMiddle(double xf, double yf, double zf,
		     double xm, double ym, double zm,
		     double xl, double yl, double zl, double Bz  );

  int createXYLast(double xf, double yf, double zf,
		   double xm, double ym, double zm,
		   double xl, double yl, double zl, double Bz  );

  int createZ( double xf, double yf, double zf,
	       double xm, double ym, double zm,
	       double xl, double yl, double zl, double Bz  );

  int getDistanceAtXY( double x, double y, double z, double &duz, double &dv, double Bz ) const;

  int getDistanceAtZ( double x, double y, double z, double &duz, double &dv, double Bz) const;
  
  double getPtkG() const { return pt;}//*Geo::OriginBzkG; }

  int getPhiZatR( double r, double &phi, double &z, double Bz ) const;
  int getPhiRatZ( double z, double &phi, double &r, double Bz ) const;

  
  void Print() const;

  static int estmateBzkG( const HitMC &mc, double x, double y, double &BzkG );

};

#endif
