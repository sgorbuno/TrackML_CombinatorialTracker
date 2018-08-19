#include "TrackModelPhysical.h"
#include <math.h>
#include <iostream>
#include "Tracker.h"


using namespace std;

int TrackModelPhysical::createXY(double xf, double yf, double zf,
				 double xm, double ym, double zm,
				 double xl, double yl, double zl, double Bz  )
{
  const double radMax = 100./Geo::OriginBzkG;//Bz; // 100 GeV maximum for pu

  x0 = xf;
  y0 = yf;
  z0 = zf;

  xm-= xf;
  ym-= yf;
  zm-= zf;

  xl-= xf;
  yl-= yf;
  zl-= zf;

  double L = sqrt( xl*xl + yl*yl );
  
  //if( L < 1. ) return -1; // 1cm cut on length
  //if( L < .1 ) return -1; // 1cm cut on length

  double Linv = 1./L;

  ex = xl*Linv; 
  ey = yl*Linv; 


  double um =  xm*ex + ym*ey;
  double vm = -xm*ey + ym*ex;

  //SG!! if( um<0.1*L || um>0.9*L ) return -2;
  //if( um<0.1 || um>L-0.1 ) return -2;

  q = 1;
  pv = 0.5*L;

  if( vm<0. ){
    pv = -pv;
    q = -1;    
  }
      
  double tmp = 0.5* ( um*(L-um) - vm*vm );
  double vmAbs = fabs(vm);
  if( tmp < 0.02*fabs(pv*vm) ) return -3; // incl. angle > 88.8 grad
  
  if( tmp > radMax*vmAbs ){
    pu = radMax; // pu > 100 GeV -> pu = 100GeV 
    pv = 0;
  } else pu = tmp/vmAbs;

  int ret = 0;

  pt = sqrt( pu*pu + pv*pv );
  //if( pt<.01 ) ret = -2;

  ptInv = 1./pt;
  pz = 0.5*zl / asin(0.5*L*ptInv);

  pu*=Bz;
  pv*=Bz;
  pz*=Bz;
  pt*=Bz;
  ptInv = 1./pt;
  return ret;
}


int TrackModelPhysical::createXYMiddle(double xf, double yf, double zf,
				       double xm, double ym, double zm,
				       double xl, double yl, double zl, double Bz  )
{
  const double radMax = 100./Geo::OriginBzkG;///Bz; // 100 GeV maximum for pu

  x0 = xm;
  y0 = ym;
  z0 = zm;

  xf-= x0;
  yf-= y0;
  zf-= z0;

  xl-= x0;
  yl-= y0;
  zl-= z0;

  double L = sqrt( xl*xl + yl*yl );
  
  //if( L < 1. ) return -1; // 1cm cut on length
  //if( L < .1 ) return -1; // 1cm cut on length

  double Linv = 1./L;

  ex = xl*Linv; 
  ey = yl*Linv; 

  double uf =  xf*ex + yf*ey;
  double vf = -xf*ey + yf*ex;

  q = 1;
  pv = 0.5*L;

  if( vf > 0. ){
    pv = -pv;
    q = -1;    
  }
      
  double tmp = 0.5* ( uf*(-L+uf) + vf*vf );
  double vfAbs = fabs(vf);
  if( tmp < 0.02*fabs(pv*vf) ) return -3; // incl. angle > 88.8 grad
  
  if( tmp > radMax*vfAbs ){
    pu = radMax; // pu > 100 GeV -> pu = 100GeV 
    pv = 0;
  } else pu = tmp/vfAbs;

   int ret = 0;

  pt = sqrt( pu*pu + pv*pv );
  //if( pt<.01 ) ret = -2;

  ptInv = 1./pt;
  pz = 0.5*zl / asin(0.5*L*ptInv);

  pu*=Bz;
  pv*=Bz;
  pz*=Bz;
  pt*=Bz;
  ptInv = 1./pt;
 return ret;
}


int TrackModelPhysical::createXYLast(double xf, double yf, double zf,
				     double xm, double ym, double zm,
				     double xl, double yl, double zl, double Bz  )
{
  const double radMax = 100./Geo::OriginBzkG;///Bz; // 100 GeV maximum for pu

  xf-= xl;
  yf-= yl;
  zf-= zl;

  xm-= xl;
  ym-= yl;
  zm-= zl;

  double L = sqrt( xm*xm + ym*ym );
  
  //if( L < 1. ) return -1; // 1cm cut on length
  //if( L < .1 ) return -1; // 1cm cut on length

  double Linv = 1./L;

  double ex1 = -xm*Linv; 
  double ey1 = -ym*Linv; 

  double uf =  xf*ex1 + yf*ey1;
  double vf = -xf*ey1 + yf*ex1;

  double tmp = 0.5* ( uf*(L+uf) + vf*vf );
  double vfAbs = fabs(vf);
  if( tmp < 0.02*fabs(pv*vf) ) return -3; // incl. angle > 88.8 grad
  
  x0 = xl;
  y0 = yl;
  z0 = zl;
  ex = ex1;
  ey = ey1;

  q = 1;
  pv = -0.5*L;

  if( vf > 0. ){
    pv = -pv;
    q = -1;    
  }
  
  if( tmp > radMax*vfAbs ){
    pu = radMax; // pu > 100 GeV -> pu = 100GeV 
    pv = 0;
  } else pu = tmp/vfAbs;
  
  pt = sqrt( pu*pu + pv*pv );
  //if( pt<.01 ) ret = -2;

  ptInv = 1./pt;
  pz = -0.5*zm / asin(0.5*L*ptInv);

  pu*=Bz;
  pv*=Bz;
  pz*=Bz;
  pt*=Bz;
  ptInv = 1./pt;
  return 0;
}
 
int TrackModelPhysical::createXYFirst(double xf, double yf, double zf,
				      double xm, double ym, double zm,
				      double xl, double yl, double zl, double Bz  )
{
  const double radMax = 100./Geo::OriginBzkG;///Bz; // 100 GeV maximum for pu

  x0 = xf;
  y0 = yf;
  z0 = zf;

  xl-= x0;
  yl-= y0;
  zl-= z0;

  xm-= x0;
  ym-= y0;
  zm-= z0;

  double L = sqrt( xm*xm + ym*ym );
  
  //if( L < 1. ) return -1; // 1cm cut on length
  //if( L < .1 ) return -1; // 1cm cut on length

  double Linv = 1./L;

  ex = xm*Linv; 
  ey = ym*Linv; 

  double ul =  xl*ex + yl*ey;
  double vl = -xl*ey + yl*ex;

  q = 1;
  pv = 0.5*L;

  if( vl > 0. ){
    pv = -pv;
    q = -1;    
  }
      
  double tmp = 0.5* ( ul*(-L+ul) + vl*vl );
  double vfAbs = fabs(vl);
  if( tmp < 0.02*fabs(pv*vl) ) return -3; // incl. angle > 88.8 grad
  
  if( tmp > radMax*vfAbs ){
    pu = radMax; // pu > 100 GeV -> pu = 100GeV 
    pv = 0;
  } else pu = tmp/vfAbs;

   int ret = 0;

  pt = sqrt( pu*pu + pv*pv );
  //if( pt<.01 ) ret = -2;

  ptInv = 1./pt;
  pz = 0.5*zm / asin(0.5*L*ptInv);

  pu*=Bz;
  pv*=Bz;
  pz*=Bz;
  pt*=Bz;
  ptInv = 1./pt;
 return ret;
}
 

int TrackModelPhysical::createZ(double xf, double yf, double zf,
				double xm, double ym, double zm,
				double xl, double yl, double zl, double Bz  )
{
  return createXYFirst(xf,yf,zf,xm,ym,zm,xl,yl,zl,Bz);


  const double puMax = 100./Geo::OriginBzkG;///Bz; // 100 GeV maximum for pu

  x0 = xf;
  y0 = yf;
  z0 = zf;

  xm-= xf;
  ym-= yf;
  zm-= zf;

  xl-= xf;
  yl-= yf;
  zl-= zf;

  double L = sqrt( xl*xl + yl*yl );
  
  if( L < 1. ) return -1; // 1cm cut on length

  double Linv = 1./L;

  ex = xl*Linv;
  ey = yl*Linv;

  double um =  xm*ex + ym*ey;
  double vm = -xm*ey + ym*ex;

  if( um<0.1*L || um>0.9*L ) return -2; //SG!! remove the cut later

  q = 1;
  pv = 0.5*L;

  if( vm<0. ){
    pv = -pv;
    q = -1;    
  }
      
  double tmp = 0.5* ( um*(L-um) - vm*vm );
  double vmAbs = fabs(vm);
  if( tmp < 0.02*fabs(pv*vm) ) return -3; // incl. angle > 88.8 grad  
  if( tmp > puMax*vmAbs ) pu = puMax; // pu > 100 GeV -> pu = 100GeV 
  else pu = tmp/vmAbs;

  int ret = 0;

  pt = sqrt( pu*pu + pv*pv );
  ptInv = 1./pt;
  //if( pt<.01 ) ret = -2;

  // calculate ptInv with respect to zm

  double cl = 0.5*L;
  double cm = 0.5*sqrt( um*um + vm*vm );
  for(int iter=0; iter<0; iter++){
    // we solve equation zm*asin(cl*ptInv) = zl*asin(cm*ptInv)
    // new ptInv = ptInv + d.
    double al = cl*ptInv;
    double am = cm*ptInv;
    // linearisation at d=0:
    // zm*( asin(al) + al/sqrt(1.-al*al)*d ) = zl*( asin(am) + am/sqrt(1.-am*am)*d )
    // d * [ zm*al/sqrt(1.-al*al) - zl*am/sqrt(1.-am*am) ] = zl*asin(am) - zm*asin(al)
    double bl = sqrt(1.-al*al);
    double bm = sqrt(1.-am*am);
    // d * [ (zm*al*bm - zl*am*bl)/(am*bm) ] = zl*asin(am) - zm*asin(al)
    double d = am*bm*( zl*asin(am) - zm*asin(al) ) / (zm*al*bm - zl*am*bl);
    ptInv+=d;
  }

  pt = 1./ptInv;
  pu = sqrt(pt*pt-pv*pv);
  pz = 0.5*zl/asin(0.5*L*ptInv);  

  pu*=Bz;
  pv*=Bz;
  pz*=Bz;
  pt*=Bz;
  ptInv = 1./pt;

  return ret;
}

int TrackModelPhysical::getDistanceAtXY( double x, double y, double z, double &duz, double &dv, double Bz ) const
{
  dv = 0;
  duz = 0;
  x-=x0;
  y-=y0;
  z-=z0;

  double u =  x*ex + y*ey;
  double v = -x*ey + y*ex;
  
  double epv = pv - q*Bz*u;
  double epu2 = pt*pt - epv*epv;
  //SG!!! if( epu2<1 ) return -1;
  
  double epu = sqrt(epu2);
  double tv = ( pv + epv )/( pu + epu );

  double vTrack = u*tv;
  
  double chord = u*sqrt(1.+tv*tv); // chord to the extrapolated point == sqrt(u^2+vTrack^2)*sign(u)
  double sa = 0.5*chord*q*Bz*ptInv; //  sin( half of the rotation angle ) ==  (chord/2) / radius

  double dS = q*pt/Bz*2.*asin( sa ); // path in XY
  double dLp = ptInv*dS; // path in XYZ / p == path in XY / pt
  double zTrack = pz * dLp;

  dv = v-vTrack;
  duz = z-zTrack;
  return 0;
}

int TrackModelPhysical::getDistanceAtZ( double x, double y, double z, double &duz, double &dv, double Bz ) const
{
  duz = 0;
  dv = 0;
  x-=x0;
  y-=y0;
  z-=z0;
 
  double u =  x*ex + y*ey;
  double v = -x*ey + y*ex;

  double angle = z/(pz/Bz); // = s/pt
  double f1 = sin(angle);
  double f2 = q*(1.-cos(angle));
  
  double uTrack = (f1*pu + f2*pv)/Bz;
  double vTrack = (f1*pv - f2*pu)/Bz;

  duz = u - uTrack;
  dv = v - vTrack;

  return 0;
}



int TrackModelPhysical::getPhiZatR( double r, double &phi, double &z, double Bz ) const
{
  // TODO: simplify

  // xc,yc center of circle

  double px = (pu*ex - pv*ey)/Bz;
  double py = (pu*ey + pv*ex)/Bz;

  double xc =  x0 + py*q;
  double yc =  y0 - px*q;

  double rc = sqrt(xc*xc+yc*yc);
  //cout<<"circle xc "<<xc<<" yc "<<yc<<" rc "<<rc<<endl;  
  if( rc<1 ) return  -1; // center of helix is too close to the origin

  double rci = 1./rc;
  double r1 = r*rci;
  double pt1 = (pt/Bz)*rci;
  double a = 0.5*( 1.  + r1*r1 - pt1*pt1 );
  
  double b = r1*r1-a*a;
  if( b<1.e-8 ) return -2;
  b = sqrt(b);

  double x = xc*a - yc*b;
  double y = yc*a + xc*b;

  double x2 = xc*a + yc*b;
  double y2 = yc*a - xc*b;
  
  if( (x-x0)*(x-x0)+(y-y0)*(y-y0) > (x2-x0)*(x2-x0)+(y2-y0)*(y2-y0) ){
    x = x2;
    y = y2;
  }
 

  double u = (x-x0)*ex + (y-y0)*ey;
  //double v =-(x-x0)*ey + (y-y0)*ex;
  
  phi = atan2(y,x);

  // get z at u:
  
  double chord = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)); // chord to the extrapolated point == sqrt(u^2+vTrack^2)*sign(u)

  if( u<0 ) chord = -chord;
  double sa = 0.5*chord*q*(ptInv*Bz); //  sin( half of the rotation angle ) ==  (chord/2) / radius
  double dS = q*(pt/Bz)*2.*asin( sa ); // path in XY
  double dLp = (ptInv*Bz)*dS; // path in XYZ / p == path in XY / pt
  z = z0 + (pz/Bz) * dLp;

  return 0;
}


int TrackModelPhysical::getPhiRatZ( double z, double &phi, double &r, double Bz ) const
{
  z-=z0;     
  int err=0;
  double angle = 0;

  if( fabs(pz/Bz)>1.e-8 ) angle = z/(pz/Bz); // = s/pt
  else err=-1;

  double f1 = sin(angle);
  double f2 = q*(1.-cos(angle));  

  double u = (f1*pu + f2*pv)/Bz;
  double v = (f1*pv - f2*pu)/Bz;
  
  double x = x0 + u*ex - v*ey;
  double y = y0 + u*ey + v*ex;
  r = sqrt(x*x+y*y);
  phi = atan2(y,x);
  return err;
}


void TrackModelPhysical::Print() const
{
  cout<<" x "<< x0<<" y "<<y0<<" z "<<z0;
  cout<<" ex "<<ex<<" ey "<<ey<<" pt "<<pt<<" pz "<<pz<<endl;
}

int TrackModelPhysical::estmateBzkG( const HitMC &mc, double x, double y, double &BzkG )
{  
  if( mc.q!=1. && mc.q!=-1. ){
    cout<<"estimate Bz: wrong charge!!!"<<endl;
    exit(1);
  }
  double ex = mc.px/mc.pt; 
  double ey = mc.py/mc.pt;
  
  x-=mc.x;
  y-=mc.y;

  double u =  x*ex + y*ey;
  double v = -x*ey + y*ex;
  double l2 = u*u + v*v;

  if( l2<.3*.3 ) return -1; // require at least 3 mm distance in XY

  BzkG = (-2*mc.q*mc.pt*v) / l2 / Geo::CLight;
  return 0;
}
