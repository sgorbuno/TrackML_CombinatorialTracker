/*
  TNtuple *nt = (TNtuple*) gDirectory->FindObjectAny("mcTracklets");  
  "part:nhits:pt:p:w:nl:fitpt:layer:x:y:z:r:phi:ephiL:ezL:ephi:ez:dv:dz:dupl");

*/

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


void TrackletConstructor::Reconstruct( const CreationCuts &crCuts, const ProlongationCuts &prCuts )
{  

  //cout<<"----------------- tracklet constructor .."<<endl; 
 
  mCreationCuts = crCuts;
  mProlongationCuts = prCuts;

  int statFitErrors=0;

  int *mcParticlesNRec3=nullptr;
  if( mTracker->mcFlag ){
    mcParticlesNRec3 = new int[mTracker->mParticles.size()];
    for( int i=0; i<mTracker->mParticles.size(); i++ ){
      mcParticlesNRec3[i] = 0;
    }
  }

  TStopwatch timer;
  
  FitLayer layers[3];
  for( int i=0; i<3; i++ ){
    FitLayer &layer = layers[i];
    layer.mTracker = mTracker;
  }
  

  LayerCreationCuts cuts1 = mCreationCuts.layerCuts[0];
  LayerCreationCuts cuts2 = mCreationCuts.layerCuts[1];
  LayerCreationCuts cuts3 = mCreationCuts.layerCuts[2];

  int iLayer1 = cuts1.iLayer;
  int iLayer2 = cuts2.iLayer;
  int iLayer3 = cuts3.iLayer;

  layers[0].Create( iLayer1, 0, 2*cuts2.cutPhiL, 2*cuts2.cutTL ); // make same bins as in layer2
  layers[1].Create( iLayer2, 0, 2*cuts2.cutPhiL, 2*cuts2.cutTL );
  layers[2].Create( iLayer3, 0, 2*cuts3.cutPhiL, 2*cuts3.cutTL );

  FitLayer &layer1 = layers[0];
  FitLayer &layer2 = layers[1];
  FitLayer &layer3 = layers[2];  
  
  int nTracklets2 = 0;
  int nTracklets2Good = 0;

  int nTracklets3 = 0;
  int nTracklets3Good = 0;

  long int statN1=0;
  long int statN2=0;
  long int statN3=0;


  // create layers for fit
  {
    for( int vol=0; vol<Geo::NVolumes; vol++){
      const LayerProlongationCuts &cuts = mProlongationCuts.volumeCuts[vol];
      Volume &v = Geo::volumes[vol];
      for( int ilayer=0; ilayer<v.nLayers; ilayer++ ){
	int id = Geo::getLayerID(vol,ilayer) ;
	FitLayer &layer = mTracker->mFitLayers[ id ];
	layer.Create( id, 0, 2*cuts.cutPhi, 2*cuts.cutT );
      }
    }
  }


  for( int ih1=0; ih1<layer1.mFitHits.size(); ih1++ ){

    statN1++;
    
    FitLayerHit h1 = layer1.mFitHits[ih1];  
    int mcID1 = -1;

    // skip cloned hits in first and last phi search bins
    if( h1.isClone ) continue;
    if( mTracker->mcFlag ) mcID1 = mTracker->mHitsMC[h1.id].partID;
        
    double extr2T;
    double extr2Phi = h1.phi;
    
    if( crCuts.useVertex ){
      if( layer2.mType==0 ){
	extr2T = layer2.mR * h1.z/h1.r;
      } else {
	extr2T = layer2.mZ * h1.r/h1.z;
      }
    } else {
      extr2T = h1.t;
    }
    
    // cut on the layer size in Z/R

    if( extr2T < layer2.mTmin - cuts2.cutTL || extr2T > layer2.mTmax + cuts2.cutTL ) continue;
    
    // look for the second hit 
    
    int iBin2Phi = layer2.getBinPhi( extr2Phi );	 
    int iBin2T = layer2.getBinT( extr2T );      
    if( iBin2Phi<0 || iBin2T<0 ) continue;

    for( int it2=0; it2<2; it2++, iBin2T++ ){
      int iBin2 = iBin2T*layer2.mNbinsPhi + iBin2Phi;
      int ih2 = layer2.mBins[iBin2].firstHit;
      int ih2end = layer2.mBins[iBin2+1].firstHit+layer2.mBins[iBin2+1].nHits;
      for( ; ih2<ih2end; ih2++){

	statN2++;

 	FitLayerHit &h2 = layer2.mFitHits[ih2];   
 
	// accurate cut on z and phi
	if( fabs( extr2T   - h2.t   ) > cuts2.cutTL   ) continue;
	if( fabs( extr2Phi - h2.searchPhi ) > cuts2.cutPhiL ) continue;

	nTracklets2++;
	
	// check mc info
	int mcID2=-1;
	if( mTracker->mcFlag ) {
	  HitMC &mc2 = mTracker->mHitsMC[h2.id];
	  mcID2 = mc2.partID;
	  if( mcID1 >=0 && mc2.partID != mcID1  ){
	    mcID2 = -1;
	  }
	  if( mcID2>=0 ) nTracklets2Good++;	  
	}

	// look for the third hit
	
	double extr3T = 0;
	double extr3Phi = 0;
	
	// make straight helix to use existing helix utility for layer crossing (TODO: use simple math )
		
	TrackModelPhysical t2Line;
	{
	  int err = t2Line.createXYLast( h1.x, h1.y, h1.z,
					 (h1.x+h2.x)/2, (h1.y+h2.y)/2,(h1.z+h2.z)/2,
					 h2.x, h2.y, h2.z, (h1.BzFwd+h2.BzBck)/2   );	    
	  if( layer3.mType == 0 ){ // radial layer
	    err = t2Line.getPhiZatR( layer3.mR, extr3Phi, extr3T, h2.BzFwd );
	  } else {
	    err = t2Line.getPhiRatZ( layer3.mZ, extr3Phi, extr3T, h2.BzFwd );
	  }
	  if( err!=0 ){
	    cout<<"Error by tracklet line extrapolation to layer "<<iLayer3<<" r: "<<layer3.mR<<" z "<<layer3.mZ<<" !!!:  err= "<<err<<endl;	    
	    layer3.Print();
	    statFitErrors++;
	    continue;
	  }
	}
	
	// cut on z   
 
	if( extr3T < layer3.mTmin - cuts3.cutTL  || extr3T > layer3.mTmax + cuts3.cutTL ) continue;
	
	int iBin3Phi = layer3.getBinPhi( extr3Phi );	 
	int iBin3T = layer3.getBinT( extr3T );      
	if( iBin3Phi<0 || iBin3T<0 ) continue;
  	for( int it3=0; it3<2; it3++, iBin3T++ ){
	  int iBin3 = iBin3T*layer3.mNbinsPhi + iBin3Phi;
	  int ih3 = layer3.mBins[iBin3].firstHit;
	  int ih3end = layer3.mBins[iBin3+1].firstHit+layer3.mBins[iBin3+1].nHits;
	  for( ; ih3<ih3end; ih3++){
	  
	    statN3++;

	    FitLayerHit &h3 = layer3.mFitHits[ih3];

 	    // accurate cut on z and phi
	    if( fabs( h3.t   - extr3T   ) > cuts3.cutTL   ) continue;
	    if( fabs( h3.searchPhi - extr3Phi ) > cuts3.cutPhiL ) continue;


	    // all cuts passed, try to create a helix
	    
	    double Bz = h2.BzMid;
	    TrackModelPhysical t;
	    
	    int err = t.createXYLast( h1.x, h1.y, h1.z,
				      h2.x, h2.y, h2.z,
				      h3.x, h3.y, h3.z, Bz);	    
	    if( err!=0 ){
	      // cout<<"Error by tracklet creation "<<" err= "<<err<<endl;
	      statFitErrors++;
	      continue;
	    }

	    // get Z deviation for the first hit
	    
	    double dv=0, duz=0;	  
	    if( layer1.mType==0  ){
	      err = t.getDistanceAtXY( h1.x, h1.y, h1.z, duz, dv, Bz);
	    } else {
	      err = t.getDistanceAtZ( h1.x, h1.y, h1.z, duz, dv, Bz);
	    }
	    if( err!=0 ){
	      cout<<"Error by track chi2 calculation!!!: err= "<<err<<endl;
	      statFitErrors++;
	      continue;
	    }
 	    
	    // cut on dz at layer 2
	    if( fabs(duz) > cuts3.cutDuz ) continue;	     
	    
	    if( t.getPtkG() < mCreationCuts.ptCut ) continue;
 
	    nTracklets3++;	    	    	    
	    int mcID3=-1;
	    if( mTracker->mcFlag ) {
	      HitMC &mc3 = mTracker->mHitsMC[h3.id];
	      if( mcID2 >=0 && mc3.partID==mcID2 ){
		mcID3 = mcID2;
		nTracklets3Good++;
		mcParticlesNRec3[mcID3]++;
	      }
	    }
	    //if( mcID3<0 ) continue;//SG!!!
	    //t.Print();
	    ProlongateTracklet( iLayer1, h1,
				iLayer2, h2,
				iLayer3, h3,
				t, vTracklets );	    
	    
	  } // ih3
	} // iz3	
      } // ih2
    } // iz2
    
  } // ih1

  timer.Stop();
 
  int nParticlesRec3=0;
  if( mTracker->mcFlag ){
    for( int i=0; i<mTracker->mParticles.size(); i++ ){
      Particle &p = mTracker->mParticles[i];
      if( p.w>0 && mcParticlesNRec3[i]>0 ) nParticlesRec3++;
    }
  }
  delete[] mcParticlesNRec3;

  cout<<"Created "<<nTracklets3<<" triplets  ("<<  nTracklets3Good <<" good )"<<endl;
  cout<<"Reconstructed "<<nParticlesRec3<<" particles at layer 3"<<endl;
  /*
  cout<<"\nfit errors: "<<statFitErrors<<endl;

  cout<<"\ntime "<<timer.CpuTime()*1000<<" ms\n"<<endl;

  cout<<"entries in 1 loop: "<<statN1<<endl;
  cout<<"entries in 2 loop: "<<statN2<<endl;
  cout<<"entries in 3 loop: "<<statN3<<endl;  
  */
}


