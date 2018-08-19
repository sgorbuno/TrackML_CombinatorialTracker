#ifndef GEO_H
#define GEO_H


#include "util.h"
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

struct Layer
{
  int id = -1; // global layer id
  int size = -1; // radial position: 0 = inner, 1 = middle, 2 = outer 
  int type = -1; // 0 = radial modules, 1 = vertical modules
  int volume = -1; //  volume number
  int layer = -1; // layer number withing the volume  
  double r = 0; 
  double z = 0;
  double tMin = 0;
  double tMax = 0;  
  
  static constexpr int NFieldPar=9;
  double field[3][NFieldPar]; // Bz field in kGaus* speed of light   (t,phi) 
  
  double getFieldFwd( double phi, double t ){ return getField( 0, phi, t ); }
  double getFieldMid( double phi, double t ){ return getField( 1, phi, t ); }
  double getFieldBck( double phi, double t ){ return getField( 2, phi, t ); }

  double getField( int idirection, double phi, double t){
    double t2 = t*t;
    double *c = field[idirection];
    return (   (c[0] + c[1]*t + c[2]*t2) + 
	       (c[3] + c[4]*t + c[5]*t2)*sin(phi) + 
	       (c[6] + c[7]*t + c[8]*t2)*cos(phi)   );
  }

};

struct Volume
{
  void init( int s, int t, int nL, int &nLayersTotal );
  int size = -1; // radial position: 0 = inner, 1 = middle, 2 = outer 
  int type = -1; // 0 = radial modules, 1 = vertical modules
  int nLayers = -1; // number of layers [2..7]
  int layerIDs[7] = {-1,-1,-1,-1,-1,-1,-1};
};

class Geo
{
 public:

  static constexpr int NVolumes = 9;
  static constexpr int NLayers = (7+4+7)+(6+4+6)+(6+2+6); // 48

  static constexpr double CLight = 0.000299792458; // speed of light 
  static constexpr double OriginBzkG = 20.*CLight; // 20 kG * speed of light 

  static Volume volumes[NVolumes];
  static Layer  layers[NLayers];

  static void init();

  static int getLayerID( int volume, int layer )  { 
    return volumes[volume].layerIDs[layer];
  }
  
  static int getVolume( int layerID )  {
    return layers[layerID].volume;
  }

  static int getVolumeLayer( int layerID )  {
    return layers[layerID].layer;
  }
 
};

#endif
