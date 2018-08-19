
#include "FitLayer.h"

#include "util.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "Tracker.h"

using namespace std;

constexpr double PI = 3.14159265359;

void  FitLayerHit::Print() const
{
  cout<<"fit hit x "<<x<<" y "<<y<<" z "<<z<<" r "<<r<<" phi "<<phi<<" t "<<t<<endl;
}

void FitLayer::Print()
{
  cout<<"FitLayer v "<<mVolume<<" l "<< mLayer<<" r "<<mR<<" z "<<mZ<<endl;
}


void FitLayer::Create( int layerID, bool phiWindowInCm, double searchWindowPhi, double searchWindowT, bool prn )
{
  //cout<<"Layer create : layer "<<layerID<<endl;
 
  // SG!!!! remove some hits !!!!

  Layer &geoLayer = Geo::layers[layerID];   
  bool *hitOK = new bool[ mTracker->mHits.size() ];  

  for( int ih=0; ih<mTracker->mHits.size(); ih++ ){    
    hitOK[ih] = 1;    
    if( !mTracker->mcFlag ) continue;
    const HitMC &mc = mTracker->mHitsMC[ih];

    //hitOK[ih] = 1;// ( mc.partID>=0 && mTracker->mParticles[mc.partID].pt>=1. );
    //hitOK[ih] = ( mc.partID>=0 && 0.5<=pt0 && pt0 <1 );
    //hitOK[ih] = ( mc.partID>=0 && pt0>=0.3 && pt0<0.5 );
    //hitOK[ih] = !( mc.partID>=0 && pt0>=0.5 );
    
    // if( mc.partID<0 ) hitOK[ih] = 0;
    
    if( mc.partID>=0 ){
      Particle &p = mTracker->mParticles[mc.partID];
      //if( p.pt >= 0.5 ) hitOK[ih] = 0;     
      //if( p.pt < 0.3 || p.pt >= 0.5 ) hitOK[ih] = 0;     
      //if( p.baseV ) hitOK[ih] = 0;
      //if( !p.baseA0 ) hitOK[ih] = 0;
      //if( p.pt>=0.3 ) hitOK[ih] = 0;
      //if( p.pt<0.2 ) hitOK[ih] = 0;
      // if( p.baseI==0 && p.pt>=0.2 ) hitOK[ih] = 0;
      //if( p.baseI==0 && p.pt<0.3 ) hitOK[ih] = 0;
      //if( !( p.pt>=0.5 && p.baseI0==0 && !p.prim) ) hitOK[ih] = 0;
      //if( p.pt<0.3 && p.baseI0==0 ) hitOK[ih] = 0;
      //if( p.baseI0!=0 ) hitOK[ih] = 0;
    }
    //if( mc.partID!=498 ) hitOK[ih] =0;
  }


  // find detector geometry from hits
  
  mType = Geo::layers[layerID].type;
  mVolume = Geo::layers[layerID].volume;
  mLayer = Geo::layers[layerID].layer;
  mTimeMark = 0;

  mTmin = mTmin = geoLayer.tMin;
  mTmax = mTmax = geoLayer.tMax;
  mR = 0;
  mZ = 0;
  double rmin = 0;
  int nHits = 0;
  for( int ih=0; ih<mTracker->mHits.size(); ih++ ){
    const Hit &h = mTracker->mHits[ih];    
    if( h.layerID!= layerID ) continue;
    if( h.isUsed!=0 ) continue;

    double t = h.t;    
    if( mTmin > t ) mTmin = t;
    if( mTmax < t ) mTmax = t;
    if( rmin > h.r ) rmin = h.r;
    mR+= h.r;
    mZ+= h.z;
    nHits++;
  }
  if( nHits>0 ){
    mR/=nHits;
    mZ/=nHits;
  } else {
    mR = 10;
    mZ = 10;
    rmin = 10;
  }

  mR = geoLayer.r;
  mZ = geoLayer.z;
  
   /*
  cout<<"layer "<<layerID<<" geo min/maxT: "<<geoLayer.tMin<<" "<<geoLayer.tMax
      <<" fitlayer: "<<mTmin<<" "<<mTmax<<endl;
  */
  if( phiWindowInCm ){
    searchWindowPhi/=rmin;    
  }

  // add some margin for T and PI, that exact max values stays inside the intervals
  mTmax+=.001;
  mPhiMin = -PI;
  double phiMax = PI+ 0.00001;
  mNbinsPhi = int( (phiMax - mPhiMin)/searchWindowPhi );
  mNbinsT = int( (mTmax - mTmin)/searchWindowT );  
  if( mNbinsPhi < 1 ) mNbinsPhi = 1;
  if( mNbinsT < 1 ) mNbinsT = 1;

  mBinPhi = (phiMax - mPhiMin)/mNbinsPhi;
  mBinT = (mTmax - mTmin)/mNbinsT;

  mBinPhiInv = 1./mBinPhi;
  mBinTInv = 1./mBinT;

  mNbinsPhi+=2;
  mNbinsT+=2;

  // 

  mNbinsTotal = mNbinsPhi*mNbinsT;

  mBins.resize(mNbinsTotal);  

  for( int ib=0; ib<mNbinsTotal; ib++ ){
    mBins[ib].firstHit=0;
    mBins[ib].nHits=0;
  }
  
  vector<int> vHitIDs;
  vector<int> vBinInd;

  for( int ih=0; ih<mTracker->mHits.size(); ih++ ){
    const Hit &h = mTracker->mHits[ih];
    if( h.layerID!= layerID ) continue;    
    if( !hitOK[ih] ) continue; //SG!!!
    if( h.isUsed!=0 ) continue;

    int it = int ( (h.t - mTmin)*mBinTInv ) ;
    int iphi = int ( (h.phi - mPhiMin )*mBinPhiInv );
         
    if( it<0 || it>=mNbinsT-2 || iphi<0 || iphi>= mNbinsPhi-2 ){
      cout<<"fit layer: wrong indexing: mark 1"<<endl;
      cout<<" it "<<it<<" of "<<mNbinsT<<" iphi "<<iphi<<" of "<<" "<<mNbinsPhi<<endl;
      cout<<"phi "<<h.phi<<" +pi="<< h.phi + PI <<endl;
      cout<<"t "<<h.t<<" tmin "<<mTmin<<" tmax "<<mTmax<<" tbin "<<mBinT<<endl;
      cout<<"it float = "<<(h.t - mTmin)*mBinTInv<<endl;
      exit(0);
    }
     
    it+=1;
    iphi+=1;
 
    if(prn) cout<<"create layer: hit at bins "<<iphi<<" "<<it<<endl;
   
    vHitIDs.push_back(ih);
    vBinInd.push_back( it*mNbinsPhi + iphi );
    
    // duplicate hits at the phi edges (phi=+-PI)

    if( h.phi > PI-mBinPhi ){ // put to the 0 phi bin
      if(prn) cout<<"create layer: add clone hit at bins "<<0<<" "<<it<<endl;
      vHitIDs.push_back(ih);
      vBinInd.push_back( it*mNbinsPhi + 0 );   
    } 
    else if( h.phi < -PI + mBinPhi ){ // put to the last phi bin
      if(prn) cout<<"create layer: add clone hit at bins "<<mNbinsPhi-1<<" "<<it<<endl;
      vHitIDs.push_back(ih);
      vBinInd.push_back( it*mNbinsPhi + mNbinsPhi-1 );   
    }    
  }
  
  nHits = vHitIDs.size();
  for( int ih=0; ih<nHits; ih++ ){
    FitLayerBin &bin = mBins[vBinInd[ih]];
    bin.nHits++;
  }
  
  int nHits1=0;
  for( int ib=0; ib<mNbinsTotal; ib++ ){
    FitLayerBin &bin = mBins[ib];
    bin.firstHit=nHits1;
    nHits1+=bin.nHits;
    bin.nHits = 0;
  }
   
  if(prn) cout<<"n hits to store:"<<nHits<<" "<<nHits1<<endl;
 
  if( nHits1!=nHits ){
    cout<<" wrong fit layer indexing: mark 2"<<endl;
    exit(0);
  }
  
  // store hits

  mFitHits.resize( nHits );

  for( int ih=0; ih<nHits; ih++ ){
    const Hit &h = mTracker->mHits[vHitIDs[ih]];
    //assert( h.layerID == layerID );
    //assert( h.isUsed==0 );
 
    FitLayerBin &bin = mBins[vBinInd[ih]];
    int ind = bin.firstHit + bin.nHits;
    bin.nHits++;
    FitLayerHit &rh = mFitHits[ind];
    rh.x = h.x;
    rh.y = h.y;
    rh.z = h.z;
    rh.r = h.r;
    rh.phi = h.phi;
    rh.searchPhi = rh.phi;
    rh.BzFwd = h.BzFwd;
    rh.BzMid = h.BzMid;
    rh.BzBck = h.BzBck;
    rh.isClone=0;
    rh.t = h.t;
    rh.id = vHitIDs[ih];
    rh.module = h.module;
    rh.timeMark=0;
  }

  // correct phi for the first and for the last phi bins
  for(int it=0; it<mNbinsT; it++){
    {
      FitLayerBin &bin = mBins[it*mNbinsPhi];
      for( int ih = 0; ih<bin.nHits; ih++){
	FitLayerHit &rh = mFitHits[bin.firstHit+ih];
	if( rh.phi < PI - mBinPhi ){
	  cout<<"wrong phi in the first search bin: phi="<<rh.phi<<endl;
	  cout<<"border "<<-PI + mBinPhi<<endl;
	  cout<<"mBinPhi "<< mBinPhi<<endl;
	  exit(1);
	}
	if( rh.phi > PI - mBinPhi){
	  rh.searchPhi = -2*PI+rh.phi;
	  rh.isClone=1;
	}
      }
    }
    {
      FitLayerBin &bin = mBins[it*mNbinsPhi + mNbinsPhi - 1];
      for( int ih = 0; ih<bin.nHits; ih++){
	FitLayerHit &rh = mFitHits[bin.firstHit + ih];
	if( rh.phi > -PI + mBinPhi  ){
	  cout<<"wrong phi in the last search bin: phi="<<rh.phi<<endl;
	  exit(1);
	}
	if( rh.phi < -PI + mBinPhi ){
	  rh.searchPhi = 2*PI+rh.phi;
	  rh.isClone=1;
	}
      }
    }
  }

  
  // check again 

  int nHits2=0;
  for( int ib=0; ib<mNbinsTotal; ib++ ){
    FitLayerBin &bin = mBins[ib];
    if( bin.firstHit!= nHits2){
      cout<<" wrong fit layer indexing: mark 3 "<<endl;
      exit(0);
    }
    nHits2+=bin.nHits;
  }
  
  if( nHits2!=nHits ){
    cout<<" wrong fit layer indexing: mark 4 "<<endl;
    exit(0);
  }
  
  delete[] hitOK;
}
