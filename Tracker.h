#ifndef TRACKER_H
#define TRACKER_H

#include "util.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "DataStructures.h"
#include "Geo.h"
#include "FitLayer.h"
#include "TrackletConstructor.h"

using namespace std;


struct CreationCuts;

struct Tracker
{
  Tracker();
  ~Tracker() = default;

  void readEvent( const char *directory, int event, bool loadMC );
  void reconstruct();
  void analyzeGeometry( bool endOfData );

  //static constexpr int NVolumes = Geo::NVolumes;
  //static constexpr int NLayers = Geo::NLayers;
  
  std::vector<FitLayer> mFitLayers = std::vector<FitLayer>(Geo::NLayers);

  std::vector<Hit> mHits;
  //  std::vector<HitRecInfo> mHitRecInfo;
  //std::vector<Track> mTracks;

  std::vector<HitMC> mHitsMC;
  std::vector<Particle> mParticles;
  TrackletConstructor proc;
  bool mcFlag;
};


#endif
