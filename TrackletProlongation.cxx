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
#include <list>

using namespace std;

  
struct Branch
{  
  int deep;
  int region;
  int layer;
  int skipHitsTimeMark;

  vector<int> vHits;
  int nLayers; 
  double chi2;
  int gapLength;
  int isStopped;
 
  TrackModelPhysical t;

  FitLayerHit h1;
  FitLayerHit h2;
  FitLayerHit h3;
  int l1;
  int l2;
  int l3;
};


struct FitRegion
{
  int lStart;
  int lEnd;
  int lIncr;
};

void TrackletConstructor::ProlongateTracklet( int layer1, FitLayerHit h01,
					      int layer2, FitLayerHit h02,
					      int layer3, FitLayerHit h03,
					      TrackModelPhysical t0,
					      vector<Tracklet> &tracklets )
{
  bool prn=0;
   
  if( layer2<0 || layer3<0 ){
    cout<<"ProlongateTracklet:: wrong input"<<endl;
    exit(1);
  }

  list<Branch> branches;
  {
    Branch br;
    br.vHits.clear();
    br.chi2 = 0;
    br.isStopped = 0;
    br.nLayers =0;
    br.deep = 0;
    br.region = 0;
    br.layer = layer3+1;
    br.skipHitsTimeMark = -1;
    br.t = t0;  
    br.h1 = h01;
    br.h2 = h02;
    br.h3 = h03;
    br.gapLength = 0;
    branches.push_back( br );
  }

  FitRegion fitRegions[3];
  fitRegions[0] = {layer3+1, Geo::NLayers, 1 }; // outwards
  fitRegions[1] = {layer3  , layer2-1, -1 }; // inbetween the base layers
  fitRegions[2] = {layer2-1, -1, -1 }; // inwards  


  int statFitErrors =0;

  for( list<Branch>::iterator iterbr=branches.begin(); iterbr!=branches.end(); iterbr++ ){ 

    Branch &br = *iterbr;

    bool regionInitialized = 1;
    for( ; br.region<3; br.region++ ){
      
      FitRegion &reg = fitRegions[br.region];

      TrackModelPhysical &t = br.t;      
      FitLayerHit &h1 = br.h1;
      FitLayerHit &h2 = br.h2;
      FitLayerHit &h3 = br.h3;
      
      if( !regionInitialized ){
	br.layer = reg.lStart;
	br.gapLength = 0;

	if( br.region==0 ){
	  h1 = h01;
	  h2 = h02;
	  h3 = h03;
	  t = t0;
	} else if( br.region==1 ){
	  h1 = h01;
	  h2 = h01;
	  h3 = h01;
	  t = t0;
	} else {
	  TrackModelPhysical tnew ;
	  int err = tnew.createXYMiddle( h01.x, h01.y, h01.z,			      
					 h02.x, h02.y, h02.z,
					 h03.x, h03.y, h03.z , h02.BzMid );
	  if( err!=0 ){
	    cout<<"Error by tracklet creation for inward prolongation !!!:  err= "<<err<<endl;
	    statFitErrors++;
	    br.skipHitsTimeMark = -1;
	    continue; // next region
	  }
	  h1 = h02; // first hit will be added at the next step
	  h2 = h03;
	  h3 = h03;  
	  t = tnew;
	}
      }
      regionInitialized = 0;
      
      for( ;  br.layer!=reg.lEnd; br.layer+=reg.lIncr ){
	
	int skipHitsTimeMark = br.skipHitsTimeMark;
	br.skipHitsTimeMark = -1;
	
	if( br.layer==layer3 || br.layer==layer2 ) prn=1;
	else prn = 0;

	if( br.gapLength>2 ){ 
	  br.isStopped = 1;
	  break; // the gap is too big, stop this direction
	}


	// propagate the track to the br.layer
	
	FitLayer &layer = mTracker->mFitLayers[ br.layer ];      
	const bool radialLayer = (layer.mType == 0);
	const LayerProlongationCuts &cut = mProlongationCuts.volumeCuts[layer.mVolume];

	double cutMult=1;
	//if( layer.mVolume==6 && layer.mLayer==1 )  cutMult=2;

	if( skipHitsTimeMark < 0 ) skipHitsTimeMark = layer.mTimeMark+1;

	int err=0;
	
	double extrBz = h3.BzFwd;
	if( br.region==1 ) extrBz = h2.BzMid;
	if( br.region==2 ) extrBz = h1.BzBck;

	
	if( ( (layer.mType==-1)&&(t.pz>=0)) || 
	    ( (layer.mType== 1)&&(t.pz<=0))    ){
	  if(prn) cout<<"layer not crossed: wrong z side"<<endl;
	  continue; // wrong z side
       }
	
	double extrT;
	double extrPhi;

	if( layer.mType == 0 ){ // radial layer	
	  err = t.getPhiZatR( layer.mR, extrPhi, extrT, extrBz );
	} else {
	  err = t.getPhiRatZ( layer.mZ, extrPhi, extrT, extrBz );
	}
	//cout<<"extrPhi = "<<extrPhi<<" extrBz = "<<extrBz<<endl;
	//layer.Print();

	if( err!=0 ) break;
      
	if( err!=0 ){
	  if(prn) cout<<"Error by tracklet extrapolation to layer r: "<<layer.mR<<" z "<<layer.mZ<<" !!!:  err= "<<err<<endl;
	  statFitErrors++;
	  continue;
	}

	// check layer crossing
	if( extrT < layer.mTmin - cut.cutT || extrT > layer.mTmax + cut.cutT ){	  
	  if(br.layer==layer1 || br.layer==layer2 || br.layer==layer3 ){
	    //SG!! cout<<"fit: track does not cross its base layer"<<endl;
	  } 
	  //if( extrT < layer.mTmin ) extrT = layer.mTmin;
	  //if( extrT > layer.mTmax ) extrT = layer.mTmax;
	  continue;
	}

	bool layerInnerCrossed = 0;
	if( layer.mVolume<4 && extrT >= layer.mTmin + cut.cutT + 2 && extrT <= layer.mTmax - cut.cutT - 2 ){
	  //if(prn) cout<<"layer inner crossed"<<endl;
	  layerInnerCrossed = 1;
	}
 
	//if(prn) cout<<"layer crossed at phi "<<extrPhi<<" t "<<extrT<<endl;
    
	int bestHit = -1;
	double bestChi2 = 1.e10;
      
	int iBinT = layer.getBinT( extrT );
	int iBinPhi = layer.getBinPhi( extrPhi );	 
         if( iBinT<0 || iBinPhi<0 ) continue;

	vector<int> hitsAround;
	prn=0;
	for( int it=0; it<2; it++, iBinT++ ){
	  int iBin = iBinT*layer.mNbinsPhi + iBinPhi;
	  int ih = layer.mBins[iBin].firstHit;
	  int ihend = layer.mBins[iBin+1].firstHit+layer.mBins[iBin+1].nHits;
	  for( ; ih<ihend; ih++){
	    FitLayerHit &h = layer.mFitHits[ih];
	    if( h.timeMark >= skipHitsTimeMark ) continue;
	    
	    double hitBz = h.BzBck;
	    if( br.region==1 ) hitBz = h.BzMid;
	    if( br.region==2 ) hitBz = h1.BzFwd;

	    //if( h.id==ihit3 ) prn=1;
	    //else prn=0;
	    if(prn) cout<<"hit in search bins: ih "<<ih<<" phi "<<h.phi<<" t "<<h.t<<endl;

	    // accurate cut on z and phi withing 4 search bins
	    if( fabs( h.t - extrT   ) > cut.cutT*cutMult   ){
	      if(prn){ cout<<"cut T not passed: "<<h.t  - extrT<<endl;
		h.Print();
		t.Print();
	      }
	      continue;
	    }
	    if( fabs( h.searchPhi - extrPhi ) > cut.cutPhi*cutMult ){
	      if(prn){
		cout<<"cut phi not passed: "<< h.phi - extrPhi<<endl;
		h.Print();
		t.Print();
	      }
	      continue;
	    }
	    
	    // get hit deviations from the trajectory

	    double duz=0;
	    double dv=0;
	  
	    if( radialLayer ){
	      err = t.getDistanceAtXY( h.x, h.y, h.z, duz, dv, (extrBz+hitBz)/2 );
	    } else {
	      err = t.getDistanceAtZ( h.x, h.y, h.z, duz, dv, (extrBz+hitBz)/2 );
	    }	  

	    if( err!=0 ){
	      cout<<"Error by track chi2 calculation!!!: err= "<<err<<endl;
	      statFitErrors++;
	      continue;
	    }
	    
	    
	    // strick cuts on deviation from the extrapolated trajectory
	  
	    //cout<<"duz "<<duz<<" dv "<<dv<<endl;
	    if( fabs( duz  ) > cut.cutUZ*cutMult ) continue;
	    if( fabs( dv   ) > cut.cutV*cutMult  ) continue;
	    
	    // accept the hit
	    hitsAround.push_back(ih);

	    double chi2 = duz*duz + dv*dv;
	    
	    if( bestChi2 > chi2 ){
	      bestChi2 = chi2;
	      bestHit = ih;
	    }
	  }
	}// bins

	//if( br.layer==layer3 ) prn=1;
	//else prn = 0;

	if(prn) cout<<"found hits in area: "<<hitsAround.size()<<endl;

	FitLayerHit hbest;
	hbest.id = -1;

	if( bestHit>=0 ) hbest = layer.mFitHits[bestHit];
	else bestChi2 = 0;
	
	// always take original hit on the base layers
	
	if(  bestHit<0 ){
	  if( br.layer == layer1 ) hbest = h01;
	  if( br.layer == layer2 ) hbest = h02;
	  if( br.layer == layer3 ) hbest = h03;
	}
	
	if( hbest.id < 0 ){
	  if( layerInnerCrossed ) br.gapLength++; // hit is missing
	  continue; // no good hit found
	}

	int timeMark = layer.mTimeMark++;
	Branch brNotUpdated = br;

	// update track parameters with the new hit
	
	err = 0;
	if( br.region==0 ){ // outwards
	  TrackModelPhysical tnew = t;
	  err = tnew.createXYLast( h2.x, h2.y, h2.z,			      
				   h3.x, h3.y, h3.z,
				   hbest.x, hbest.y,  hbest.z  , h3.BzMid );
	  if(err==0){
	    h1 = h2;
	    h2 = h3;
	    h3 = hbest;
	    t = tnew;
	  }
	} else if( br.region==2  ){ // inwards 
	  TrackModelPhysical tnew = t;
	  err = tnew.createXYFirst( hbest.x,  hbest.y,  hbest.z,			      
				    h1.x, h1.y, h1.z,
				    h2.x, h2.y, h2.z , h1.BzMid );
	  if( err==0 ){	    
	    h3 = h2;
	    h2 = h1;
	    h1 = hbest;
	    t = tnew;	  
	  }	  
	}
	
	if( err!=0 ){  
	  statFitErrors++;
	  continue;
	}      
	
	// store the hit
      
	br.gapLength = 0;
	br.nLayers++;
	br.chi2 += bestChi2;
	br.vHits.push_back(hbest.id);
	hbest.timeMark = timeMark;

	// pick up duplications of the new hit from the neighbouring modules
	int nRest=0;
	for( int ih=0; ih<hitsAround.size(); ih++ ){

	  FitLayerHit &h = layer.mFitHits[hitsAround[ih]];

	  if( h.id == hbest.id  ) continue; // the hit is already stored
	  
	  // get hit deviations from the updated trajectory

	  double duz=0;
	  double dv=0;
	  
	  if( radialLayer ){
	    err = t.getDistanceAtXY( h.x, h.y, h.z, duz, dv, hbest.BzMid );
	  } else {
	    err = t.getDistanceAtZ( h.x, h.y, h.z, duz, dv, hbest.BzMid );
	  }	  

	  if( err!=0 ){
	    cout<<"Error by track chi2 calculation!!!: err= "<<err<<endl;
	    statFitErrors++;
	    continue;
	  }
	  
	  if( fabs( duz  ) > cut.cutDuplUZ ||
	      fabs( dv   ) > cut.cutDuplV  ){
	    nRest++;
	    continue;
	  }
	  // store the hit
	  h.timeMark = timeMark;
	  br.vHits.push_back(h.id);
	}
	
	// store the other branch
	if( 0 && br.deep<5 && nRest>0 && bestHit >=0 ){
	  brNotUpdated.skipHitsTimeMark = timeMark;
	  brNotUpdated.deep = br.deep +1;
	  branches.push_back(brNotUpdated);
	}
	
	// go to the next layer     
	
      } // br.layers

    } // fit regions
    
    // store  the branch
    Tracklet tracklet;
    tracklet.vHits = br.vHits;
    tracklet.nLayers = br.nLayers; 
    tracklet.chi2 = br.chi2;
    tracklet.isStopped = br.isStopped;
    tracklet.pt = t0.getPtkG();
    tracklet.vtxTrack = (layer1<0);
    if( tracklet.nLayers<2 ) continue;
    tracklets.push_back(tracklet);	

  } // branches

  //cout<<" total n branches "<<branches.size()<<endl;

}

