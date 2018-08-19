
#include "TrackletConstructor.h"
#include "Tracker.h"
#include "util.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "FitLayer.h"
#include "TStopwatch.h"

using namespace std;


void TrackletConstructor::ReconstructPrim( const CreationPrimCuts &crCuts, const ProlongationCuts &prCuts )
{  
  
  cout<<"----------------- tracklet constructor .."<<endl; 
 
  mCreationPrimCuts = crCuts;
  mProlongationCuts = prCuts; 

  int statFitErrors=0;

  int *mcParticlesNRec2=nullptr;
  if( mTracker->mcFlag ){
    mcParticlesNRec2 = new int[mTracker->mParticles.size()];
    for( int i=0; i<mTracker->mParticles.size(); i++ ){
      mcParticlesNRec2[i] = 0;
    }
  }
  
  TStopwatch timer;
  
  FitLayer layers[2];
  for( int i=0; i<2; i++ ){
    FitLayer &layer = layers[i];
    layer.mTracker = mTracker;
  }
  

  LayerCreationCuts cuts1 = crCuts.layerCuts[0];
  LayerCreationCuts cuts2 = crCuts.layerCuts[1];

  int iLayer1 = cuts1.iLayer;
  int iLayer2 = cuts2.iLayer;


  layers[0].Create( iLayer1, 0, 2*cuts2.cutPhiL, 2*cuts2.cutTL ); // make bins as in layer2
  layers[1].Create( iLayer2, 0, 2*cuts2.cutPhiL, 2*cuts2.cutTL );

  FitLayer &layer1 = layers[0];
  FitLayer &layer2 = layers[1];
  
  int nTracklets2 = 0;
  int nTracklets2Good = 0;

  long int statN1=0;
  long int statN2=0;

  // create layers for fit
  
  for( int vol=0; vol<Geo::NVolumes; vol++){
    const LayerProlongationCuts &cuts = mProlongationCuts.volumeCuts[vol];
    Volume &v = Geo::volumes[vol];
    for( int ilayer=0; ilayer<v.nLayers; ilayer++ ){
      int id = Geo::getLayerID(vol,ilayer) ;
      FitLayer &layer = mTracker->mFitLayers[ id ];
      layer.Create( id, 0, 2*cuts.cutPhi, 2*cuts.cutT );
    }
  }
  
  //

  for( int ih1=0; ih1<layer1.mFitHits.size(); ih1++ ){

    statN1++;
    
    FitLayerHit h1 = layer1.mFitHits[ih1];   

    // skip cloned hits in first and last phi search bins
    if( h1.isClone ) continue;
    int mcID1 = -1;
    if( mTracker->mcFlag ) mcID1 = mTracker->mHitsMC[h1.id].partID;  	    
     
    double extr2Tslope;
    double extr2T;
    double extr2Phi =  h1.phi;
    if( layer2.mType==0 ){
      extr2Tslope = h1.z/h1.r; 
      extr2T = layer2.mR * extr2Tslope;
    } else {
      extr2Tslope = h1.r/h1.z;
      extr2T = layer2.mZ * extr2Tslope;
    }
    
    // cut on the layer size in Z/R

    if( extr2T < layer2.mTmin - cuts2.cutTL || extr2T > layer2.mTmax + cuts2.cutTL ) continue;

    // look for the second hit 
  
    int iBin2z = layer2.getBinT( extr2T );
    int iBin2phi = layer2.getBinPhi( extr2Phi );
    if( iBin2z<0 || iBin2phi<0 ) continue;
    int iBin2 = iBin2z*layer2.mNbinsPhi + iBin2phi;
 
    for( int iz2=0; iz2<2; iz2++ ){
      int ih2 = layer2.mBins[iBin2].firstHit;
      int ih2end = layer2.mBins[iBin2+1].firstHit+layer2.mBins[iBin2+1].nHits;
      for( ; ih2<ih2end; ih2++){
	statN2++;

 	FitLayerHit &h2 = layer2.mFitHits[ih2];   
 
	// accurate cut on z and phi
	if( fabs( extr2T   - h2.t   ) > cuts2.cutTL   ) continue;
	if( fabs( extr2Phi - h2.searchPhi ) > cuts2.cutPhiL ) continue;

	// more accurate cut on t
	
	double extr2THit;
	if( layer2.mType==0 ){
	  extr2THit = h2.r * extr2Tslope;
	  if( fabs( extr2THit   - h2.z   ) > cuts2.cutTH  ) continue;
 	} else {
	  extr2THit = h2.z * extr2Tslope;
	  if( fabs( extr2THit   - h2.r   ) > cuts2.cutTH  ) continue;
 	}
   
 	// all angular cuts passed, try to create a helix
	    	    
	TrackModelPhysical t;
	double Bz = h1.BzMid;
	
	int err = t.createXYLast( 0., 0., 0.,
				  h1.x, h1.y, h1.z,
				  h2.x, h2.y, h2.z , Bz  );	    
		
	if( err!=0 ){
	  cout<<"Error by tracklet creation "<<" err= "<<err<<endl;
	  statFitErrors++;
	  continue;
	}

	// get Z deviation for the vertex
	    
	double dv=0, duz=0;	  	
	err = t.getDistanceAtXY( 0, 0, 0, duz, dv, Bz );

	if( err!=0 ){
	  cout<<"Error by track chi2 calculation!!!: err= "<<err<<endl;
	  statFitErrors++;
	  continue;
	}
 	    
	// cut on dz at layer 2
	if( fabs(duz) > cuts2.cutDuz ) continue;
	//if( fabs(dv) > cuts3.cutDv ) continue;
 
	// cut on Pt
	//cout<<t.getPtkG()<<endl;
	if( t.getPtkG() < crCuts.ptCut ) continue;

	nTracklets2++;	    	    
	
	// check mc info
	int mcID2=-1;
	if( mTracker->mcFlag ){
	  HitMC &mc2 = mTracker->mHitsMC[h2.id];
	  mcID2 = -1;
	  if( mcID1 >=0 && mc2.partID == mcID1  ){
	    mcID2 = mcID1;
	    nTracklets2Good++;
	    mcParticlesNRec2[mcID2]++;
	  }
	}

	//if( mcID2<0 ) continue; //!!!SG!!!

	FitLayerHit h0;
	h0.x = 0;
	h0.y = 0;
	h0.z = 0;
	h0.r = 0;
	h0.phi=0;
	h0.t = 0;
	h0.BzFwd = Geo::OriginBzkG;
	h0.BzMid = Geo::OriginBzkG;
	h0.BzBck = Geo::OriginBzkG;

	ProlongateTracklet( -1,      h0,
			    iLayer1, h1,
			    iLayer2, h2,
			    t, vTracklets );	    	
      } // ih2
      iBin2 += layer2.mNbinsPhi;
    } // iz2
    
  } // ih1

  timer.Stop();
 
  int nParticlesRec2=0;

  if( mTracker->mcFlag ){
    for( int i=0; i<mTracker->mParticles.size(); i++ ){
      Particle &p = mTracker->mParticles[i];
      if( p.w>0 && mcParticlesNRec2[i]>0 ) nParticlesRec2++;
    }
  }
  delete[] mcParticlesNRec2;
  
  cout<<"Created "<<nTracklets2<<" duplets  ("<<  nTracklets2Good <<" good )"<<endl;
  cout<<"Reconstructed "<<nParticlesRec2<<" particles at layer 2"<<endl;
  
  /*
  cout<<"\nfit errors: "<<statFitErrors<<endl;

  cout<<"\ntime "<<timer.CpuTime()*1000<<" ms\n"<<endl;

  cout<<"entries in 1 loop: "<<statN1<<endl;
  cout<<"entries in 2 loop: "<<statN2<<endl;
  */  
}


