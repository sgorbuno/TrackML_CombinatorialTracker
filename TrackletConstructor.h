#ifndef TRACKLETCONSTRUCTOR_H
#define TRACKLETCONSTRUCTOR_H

#include "util.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include <iostream>
#include "TrackModelPhysical.h"

struct Tracker;
struct FitLayerHit;

#include "Cuts.h"

using namespace std;

struct Tracklet
{
  vector<int> vHits;
  int nLayers; 
  double chi2;
  int isStopped;
  double pt;
  bool vtxTrack;
};



struct TrackletConstructor
{
  
  void ReconstructMC( const CreationCuts &crCuts, const ProlongationCuts &prCuts );
  void Reconstruct( const CreationCuts &crCuts, const ProlongationCuts &prCuts );

  void ReconstructPrim  ( const CreationPrimCuts &crCuts, const ProlongationCuts &prCuts );
  void ReconstructPrimMC( const CreationPrimCuts &crCuts, const ProlongationCuts &prCuts );

  void ProlongateTracklet( int layer1, FitLayerHit h01,
			   int layer2, FitLayerHit h02,
			   int layer3, FitLayerHit h03,
			   TrackModelPhysical t0,
			   vector<Tracklet> &tracklets );

  void CleanTracklets( int minTrackNLayers );


  void TrackletFitTest( const CreationCuts &crCuts );
  void TrackletFitTest( const CreationPrimCuts &crCuts );
  void TrackletFitTest( int layer1, int layer2, int layer3 );


  void TrackletEfficiency( int useBaseLayers, const double ptCutMin, const double ptCutMax );

  Tracker *mTracker;  
  CreationCuts mCreationCuts;
  CreationPrimCuts mCreationPrimCuts;
  ProlongationCuts mProlongationCuts;

  vector<Tracklet> vTracklets;
  vector<Tracklet> vTracks;
};

#endif
