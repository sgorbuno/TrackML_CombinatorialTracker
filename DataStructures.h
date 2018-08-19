#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include <iostream>
#include <vector>


using namespace std;

 
struct Hit // container for hit information
{
  double x;
  double y;
  double z;
  double r;
  double phi;
  double t;
  double BzFwd; // magnetic field for fwd prolongation
  double BzMid;// magnetic field for helix creation
  double BzBck;// magnetic field for backward prolongation
  int volume; // volume number = 0..8
  int layer; // layer number withing the volume = 0..6  
  int layerID;
  int module; // module withing the layer
  int isUsed;
  int trackID;
  void Print() const;
};


struct HitMC // hit truth data, stored in an array parallel to hits
{  
  int hitID; // hit id
  double x;
  double y;
  double z;
  int partID;
  double w; // weight 
  double px;
  double py;
  double pz;
  double pt;
  double dOrigin2; // squared distance to particle origin 
  double p; // momentum for sorting
  double q;
  double r;
  double phi;
  double t;
  void Print() const;

  bool operator< ( const HitMC &h ) const
  {
    return ( p > h.p ) || (p==h.p && dOrigin2 < h.dOrigin2);
  }
};


struct Particle // structure for truth particle info
{
  Particle(int nhits=0)
    : hits(nhits)
  {
    hits.clear();
  }
  ~Particle() = default;
  void Print() const;

  double x;
  double y;
  double z;
  double r;  
  double px;
  double py;
  double pz;
  double pt;
  double p;
  double q;
  double w;
  double xl;
  double yl;
  double zl;
  double rl;

  int nLayers;
  std::vector<int> hits;
  std::vector<int> hitClusterIds; // along the trajectory
  int prim; // is it coming from the origin
  bool baseV;  // some combinations of layers it crosses
  bool baseV0; 
  bool baseV1; 
  bool baseV2; 
  bool baseA0;
  bool baseA1;
  bool baseA2;
  bool baseA3;
};


#endif
