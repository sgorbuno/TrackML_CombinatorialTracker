#ifndef FITLAYER_H
#define FITLAYER_H

#include "util.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include <iostream>
#include "TrackModelPhysical.h"

struct Tracker;

using namespace std;

struct FitLayerHit{
  double x;
  double y;
  double z;
  double r;
  double phi;
  double searchPhi;
  double t;
  double BzFwd;
  double BzMid;
  double BzBck;

  bool isClone;
  int id;
  int module;
  int timeMark;
  void Print() const;
};

struct FitLayerBin{
  int firstHit;
  int nHits;
};

struct FitLayer
{
  // phi + t (== z or r) coordinates
  
  void Create( int layerID, bool phiWindowInCm, double searchWindowPhi, double searchWindowT, bool prn=0  );  
  
  int getBinPhi( double phi ){ // get bin for phi with searchWindowPhi/2 marging
    int bin = int ( 0.5 + (phi - mPhiMin) * mBinPhiInv );
    if( bin < 0 || bin >= mNbinsPhi-1 ){
      cout<<"Wrong bin Phi : phi "<<phi<<" bin "<< bin <<" out of "<< mNbinsPhi <<endl;
      return -1;//exit(0);
    }
    return bin;
  }
  
  int getBinT( double t ){// get bin for t with searchWindowT/2 marging
    int bin = int ( 0.5 + (t-mTmin) * mBinTInv );
    // TODO: can be removed later

    if( bin < 0 ) bin=0;
    if( bin > mNbinsT-2 ) bin = mNbinsT-2;

    if( bin < 0 || bin > mNbinsT-2 ){
      cout<<"Wrong bin T "<< bin <<" out of "<< mNbinsT <<endl;
      return -1;//exit(0);
    }
    return bin;
  }

  int getBin( double phi, double t ){
    int it = getBinT(t);
    int iphi = getBinPhi( phi );
    return it*mNbinsPhi + iphi;
  }

  void Print();

  Tracker *mTracker = 0;
  int mType;
  int mVolume;
  int mLayer;
  double mBinPhi;
  double mBinT;
  double mBinPhiInv;
  double mBinTInv;
  int mNbinsPhi;
  int mNbinsT;
  int mNbinsTotal;
  double mTmin;
  double mTmax;
  double mPhiMin;

  double mR;
  double mZ;
  
  int mTimeMark;

  vector<FitLayerBin> mBins;
  vector<FitLayerHit> mFitHits;
  
};

#endif
