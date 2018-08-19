#ifndef CUTS_H 
#define CUTS_H 

#include "Geo.h"

struct LayerProlongationCuts
{
  double cutPhi;
  double cutT;
  double cutV;
  double cutUZ;
  double cutDuplV;
  double cutDuplUZ;
};

struct ProlongationCuts
{
  LayerProlongationCuts volumeCuts[Geo::NVolumes];
};

struct LayerCreationCuts
{
  int iLayer;
  double cutPhiL;
  double cutTL;
  double cutTH;
  double cutDv;
  double cutDuz;
};

struct CreationCuts
{
  int useVertex;
  double ptCut;
  LayerCreationCuts layerCuts[3];
};

struct CreationPrimCuts
{
  double ptCut;
  LayerCreationCuts layerCuts[2];
};


struct Cuts
{
  static void init();

  static CreationCuts mirror( const CreationCuts  &fwd ) ;

  // == primary
 
  static CreationPrimCuts CrCutsPrimV0Pt05;
  static ProlongationCuts PrCutsPrimV0Pt05;
  static CreationPrimCuts CrCutsPrimV0Pt03;
  static ProlongationCuts PrCutsPrimV0Pt03;
  static CreationPrimCuts CrCutsPrimV0Pt02;
  static ProlongationCuts PrCutsPrimV0Pt02;
  static CreationPrimCuts CrCutsPrimV0Pt01;
  static ProlongationCuts PrCutsPrimV0Pt01;

  static CreationPrimCuts CrCutsPrimV1Pt05;
  static CreationPrimCuts CrCutsPrimV2Pt05;
  static ProlongationCuts PrCutsPrimV1Pt05;

  static CreationPrimCuts CrCutsPrimV1Pt03;
  static CreationPrimCuts CrCutsPrimV2Pt03;
  static ProlongationCuts PrCutsPrimV1Pt03;

  static CreationPrimCuts CrCutsPrimV1Pt02;
  static CreationPrimCuts CrCutsPrimV2Pt02;
  static ProlongationCuts PrCutsPrimV1Pt02;
 
  static CreationPrimCuts CrCutsPrimV1Pt01;
  static CreationPrimCuts CrCutsPrimV2Pt01;
  static ProlongationCuts PrCutsPrimV1Pt01;

  // == secondary

  static CreationCuts     CrCutsA0Pt05;
  static ProlongationCuts PrCutsA0Pt05;
  static CreationCuts     CrCutsA0Pt03;
  static ProlongationCuts PrCutsA0Pt03;
  static CreationCuts     CrCutsA0Pt02;
  static ProlongationCuts PrCutsA0Pt02;
  static CreationCuts     CrCutsA0Pt01;
  static ProlongationCuts PrCutsA0Pt01;

  static CreationCuts     CrCutsB0Pt05;
  static ProlongationCuts PrCutsB0Pt05;
  static CreationCuts     CrCutsB0Pt03;
  static ProlongationCuts PrCutsB0Pt03;
  static CreationCuts     CrCutsB0Pt01;
  static ProlongationCuts PrCutsB0Pt01;

  static CreationCuts     CrCutsA2Pt01;
  static ProlongationCuts PrCutsA2Pt01;
  static CreationCuts     CrCutsB2Pt01;
  static ProlongationCuts PrCutsB2Pt01;
  static CreationCuts     CrCutsC2Pt01;
  static ProlongationCuts PrCutsC2Pt01;
  static CreationCuts     CrCutsD2Pt01;
  static ProlongationCuts PrCutsD2Pt01;
  static CreationCuts     CrCutsE2Pt01;
  static ProlongationCuts PrCutsE2Pt01;

  static CreationCuts     CrCutsA3Pt03;
  static ProlongationCuts PrCutsA3Pt03;
  static CreationCuts     CrCutsB3Pt03;
  static ProlongationCuts PrCutsB3Pt03;
  static CreationCuts     CrCutsC3Pt03;
  static ProlongationCuts PrCutsC3Pt03;
  static CreationCuts     CrCutsD3Pt05;
  static ProlongationCuts PrCutsD3Pt05;
  static CreationCuts     CrCutsD3Pt03;
  static ProlongationCuts PrCutsD3Pt03;
  static CreationCuts     CrCutsD3Pt01;
  static ProlongationCuts PrCutsD3Pt01;


  static CreationCuts     CrCutsA5Pt03;
  static ProlongationCuts PrCutsA5Pt03;


  //=== old cuts, some still in use


  //==========   Radial cuts
  
  static CreationCuts CrCutsPt1Rad;
  static CreationCuts CrCutsPt05Rad;

  static CreationCuts CrCutsPt05RadMiddle;
  static CreationCuts CrCutsPt05RadMiddleL;
  static CreationCuts CrCutsPt05MiddleI2;
  static CreationCuts CrCutsPt03Rad;
  static CreationCuts CrCutsPt01Rad;
  static CreationCuts CrCutsPt01MiddleRad;
  static CreationCuts CrCutsPtOld11MiddleRad;
  static CreationCuts CrCutsPt01RadMv1;
 
 
  //==========  Forward cuts

  static CreationCuts CrCutsPt1Fwd; 
  static CreationCuts CrCutsPt05Fwd;
  static CreationCuts CrCutsPt05FwdL;
  static CreationCuts CrCutsPt05FwdC;
  static CreationCuts CrCutsPt05FwdMiddle;
  static CreationCuts CrCutsPt05FwdMiddleL;
  static CreationCuts CrCutsPt05FwdMiddleI1;
  static CreationCuts CrCutsPt03Fwd;
  static CreationCuts CrCutsPt01Fwd;
  static CreationCuts CrCutsPt1MiddleFwd;
  static CreationCuts CrCutsPt01MiddleFwd;

  // ===========   prolongation cuts ====

  
  static ProlongationCuts PrCutsPt1Fwd ;
  static ProlongationCuts PrCutsPt05Fwd ;
  static ProlongationCuts PrCutsPt05FwdMiddle ;
  static ProlongationCuts PrCutsPt05FwdMiddleL ;
  static ProlongationCuts PrCutsPt05FwdMiddleI1 ;
  static ProlongationCuts PrCutsPt01Fwd ;

  static ProlongationCuts PrCutsPt05RadLoose ;
  static ProlongationCuts PrCutsPt05Rad ;
  static ProlongationCuts PrCutsPt05RadMiddle ;
  static ProlongationCuts PrCutsPt05I1;
  static ProlongationCuts PrCutsPt03Rad;
  static ProlongationCuts PrCutsPt01Rad ; 
  static ProlongationCuts PrCutsPt01MiddleRad ;
 
};


#endif
