/*
  TNtuple *seeds = (TNtuple*) gDirectory->FindObjectAny("mcTracklets");  
  
  "part:nhits:pt:p:w" + ":nl:prim:layer:dupl:x" + 
  ":y:z:r:phi:t" + ":dPL:dTL:dTH:duz:dv" +
  ":fitpt:bV:bV0:bV1:bV2:bA0:bA1:bA2:bA3";

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
#include "TMath.h"

using namespace std;


void TrackletConstructor::ReconstructMC( const CreationCuts &crCuts, const ProlongationCuts &prCuts )
{  
  mCreationCuts = crCuts;
  mProlongationCuts = prCuts;

  // loop over particles
  cout<<"\n-----------------\n tracklet constructor MC  .."<<endl;
  
  static int nEvents = 0;
  static TFile *statFile = 0;
  static TNtuple *ntMCTracklets = 0;

  if( nEvents == 0 ){
    statFile = new TFile("statMCTracklets.root", "RECREATE");  
    string vars;
    vars = vars + 
      "part:nhits:pt:p:w" + ":nl:prim:layer:dupl:x" + 
      ":y:z:r:phi:t" + ":dPL:dTL:dTH:duz:dv" +
      ":fitpt:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";
    
    ntMCTracklets = new TNtuple
      ("mcTracklets","mc tracklet layers", vars.data() );
    ntMCTracklets->SetMarkerStyle(8);
    ntMCTracklets->SetMarkerSize(0.3);
  }
  nEvents++;

  static int statCreateTries=0;
  static int statCreateErrors=0;

  int layerIDs[3];
  FitLayer layers[3];
  for( int i=0; i<3; i++ ){
    layerIDs[i] = mCreationCuts.layerCuts[i].iLayer;
    FitLayer &layer = layers[i];
    layer.mTracker = mTracker;
    if( layerIDs[i]>=0 ) layer.Create( layerIDs[i], 0, 0.2, 4. );
  } 
  
 
  
  for( unsigned int ipart=0; ipart<mTracker->mParticles.size(); ipart++ ){      

    Particle &p = mTracker->mParticles[ipart];     
    
    // check the number of crossed layers
    int partNlayers = 0;
    {
      int present[3] = {0,0,0};
      for( int ih = 0; ih<p.hits.size(); ih++ ){
	Hit &h = mTracker->mHits[p.hits[ih]];
	for( int il=0; il<3; il++ ){
	  if( h.layerID == layerIDs[il]) present[il] = 1;
	}
      } 
      if( !present[0] ) partNlayers = 0;
      else if( !present[1] ) partNlayers = 1;
      else if( !present[2] ) partNlayers = 2;
      else partNlayers = 3;
    }
    
    if( partNlayers < 3 ) continue;
    /*
      "part:nhits:pt:p:w" + ":nl:prim:layer:dupl:x" + 
      ":y:z:r:phi:t" + ":dPL:dTL:dTH:duz:dv" +
      ":fitpt:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";
    */
    float f[100];
    for( int j=0; j<100; j++ ) f[j] = 0.;
    f[ 0] = ipart;
    f[ 1] = p.hits.size();
    f[ 2] = p.pt;
    f[ 3] = p.p;
    f[ 4] = p.w;    
    f[ 5] = p.nLayers;
    f[ 6] = p.prim;

    f[21] = p.baseV;
    f[22] = p.baseV0;
    f[23] = p.baseV1;
    f[24] = p.baseV2;
    f[25] = p.baseA0;
    f[26] = p.baseA1;
    f[27] = p.baseA2;
    f[28] = p.baseA3;

    int duplicate[3] = {-1,-1,-1};

    duplicate[0] = -1;
    
    int h1N = p.hits.size();
    
    for( int ih1 = 0; ih1<h1N; ih1++ ){
      Hit h1 = mTracker->mHits[p.hits[ih1]];
      if( h1.layerID != layerIDs[0] ) continue;       
      duplicate[0]++;      
      /*
	"part:nhits:pt:p:w" + ":nl:prim:layer:dupl:x" + 
	":y:z:r:phi:t" + ":dPL:dTL:dTH:duz:dv" +
	":fitpt:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";
      */
      f[ 7] = 1; // layer
      f[ 8] = duplicate[0];
      f[ 9] = h1.x;
      f[10] = h1.y;
      f[11] = h1.z;
      f[12] = h1.r;
      f[13] = h1.phi;
      f[14] = h1.t;      
      f[15] = 0;
      f[16] = 0;
      f[17] = 0;
      f[18] = 0;
      f[19] = 0;      
      f[20] = 0;      
      
      ntMCTracklets->Fill(f);
      
      duplicate[1] = -1;
      for( int ih2 = 0; ih2<p.hits.size(); ih2++ ){
	Hit &h2 = mTracker->mHits[p.hits[ih2]];
	if( h2.layerID != layerIDs[1] ) continue;		
	duplicate[1]++;	
	/*
	  "part:nhits:pt:p:w" + ":nl:prim:layer:dupl:x" + 
	  ":y:z:r:phi:t" + ":dPL:dTL:dTH:duz:dv" +
	  ":fitpt:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";
	*/
	f[ 7] = 2; // layer
	f[ 8] = duplicate[1];
	f[ 9] = h2.x;
	f[10] = h2.y;
	f[11] = h2.z;
	f[12] = h2.r;
	f[13] = h2.phi;
	f[14] = h2.t;
	f[15] = h1.phi-h2.phi; // extrapolated to layer phi
	if( f[15] > TMath::Pi() )  f[15] -= TMath::TwoPi();
	if( f[15] < -TMath::Pi() )  f[15] += TMath::TwoPi();
	if( crCuts.useVertex ){
	  if( layers[1].mType==0 ){
	    double extrTL = layers[1].mR * h1.z/h1.r;
	    double extrTH =        h2.r  * h1.z/h1.r;
	    f[16] = h2.t-extrTL; // t diff at layer radius
	    f[17] = h2.z-extrTH; // z diff t hit radius
	  } else {
	    double extrTL = layers[1].mZ * h1.r/h1.z;
	    double extrTH =        h2.z  * h1.r/h1.z;
	    f[16] = h2.t - extrTL;// t diff at layer z	
	    f[17] = h2.r - extrTH;// r diff at hit z
	  }    
	}
	else {
	  if( layers[1].mType==0 ){
	    f[16] = h2.t-h1.t; // t diff at layer radius
	    f[17] = h2.z-h1.z; // z diff t hit radius
	  } else {
	    f[16] = h2.t - h1.t;// t diff at layer z	
	    f[17] = h2.r - h1.r;// r diff at hit z
	  }    
	}
	f[18] = 0;
	f[19] = 0;
	f[20] = 0;
	ntMCTracklets->Fill(f);
	
	double extr3Tslope = 0;
	if( layers[2].mType==0 ){	  
	  extr3Tslope =  (h2.z-h1.z)/(h2.r-h1.r);
	} else {
	  extr3Tslope =  (h2.r-h1.r)/(h2.z-h1.z);
	}
	
	double linePhiL=0, lineTL=0;
	TrackModelPhysical t2Line;
	{
	  int err = t2Line.createXYLast( h1.x, h1.y, h1.z,
					 (h1.x+h2.x)/2, (h1.y+h2.y)/2,(h1.z+h2.z)/2,
					 h2.x, h2.y, h2.z, (h1.BzFwd + h2.BzBck)/2   );	    
	  if( layers[2].mType == 0 ){ // radial layer	
	    err = t2Line.getPhiZatR( layers[2].mR, linePhiL, lineTL, h2.BzFwd  );
	  } else {
	    err = t2Line.getPhiRatZ( layers[2].mZ, linePhiL, lineTL, h2.BzFwd  );
	  }
	  if( err!=0 ){
	    cout<<"Error by tracklet line extrapolation to layer r: "<<layers[2].mR<<" z "<<layers[2].mZ<<" !!!:  err= "<<err<<endl;
	    statCreateErrors++;
	    continue;
	  }
	}      
	
	duplicate[2] = -1;
	for( int ih3 = 0; ih3<p.hits.size(); ih3++ ){
	  Hit &h3 = mTracker->mHits[p.hits[ih3]];
	  if( h3.layerID != layerIDs[2] ) continue;
	  duplicate[2]++;
	  /*
	    "part:nhits:pt:p:w" + ":nl:prim:layer:dupl:x" + 
	    ":y:z:r:phi:t" + ":dPL:dTL:dTH:duz:dv" +
	    ":fitpt:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";
	  */
	  f[ 7] = 3; // layer
	  f[ 8] = duplicate[2];
	  f[ 9] = h3.x;
	  f[10] = h3.y;	
	  f[11] = h3.z;
	  f[12] = h3.r;
	  f[13] = h3.phi;
	  f[14] = h3.t;
	  f[15] = h2.phi - linePhiL; // extrapolated to layer phi
	  if( f[15] > TMath::Pi() )  f[15] -= TMath::TwoPi();
	  if( f[15] < -TMath::Pi() )  f[15] += TMath::TwoPi();
	  double extrTL = 0, extrTH = 0;
	  if( layers[2].mType==0 ){
	    extrTL = lineTL;//h1.z + extr3Tslope*(layers[2].mR - h1.r);
	    extrTH = h1.z + extr3Tslope*(       h3.r  - h1.r);
	    f[16] = h3.t - extrTL;// extrapolated to layer T
	    f[17] = h3.z - extrTH;// extrapolated to hit T
	  } else {
	    extrTL = lineTL;//h1.r + extr3Tslope*(layers[2].mZ - h1.z);
	    extrTH = h1.r + extr3Tslope*(       h3.z  - h1.z);
	    f[16] = h3.t - extrTL;// extrapolated to layer T
	    f[17] = h3.r - extrTH;// extrapolated to hit T
	  }

	  // create tracklet	  
	  statCreateTries++;
	  TrackModelPhysical t;
	  int err = t.createXYLast( h1.x, h1.y, h1.z,			      
				    h2.x, h2.y, h2.z,
				    h3.x, h3.y, h3.z , h2.BzMid );	  
	  if( err!=0 ){
	    statCreateErrors++;
	    cout<<"Error by tracklet creation!!!: nhits="<<p.hits.size()<<" err= "<<err<<endl;
	    cout<<h1.x<<" "<<h1.y<<" "<<h1.z<<" "<<endl;
	    cout<<h2.x<<" "<<h2.y<<" "<<h2.z<<" "<<endl;
	    cout<<h3.x<<" "<<h3.y<<" "<<h3.z<<" "<<endl;
	    continue;
	  }

	  double dv=0, duz=0;	  
	  if( layers[1].mType==0  ){
	    err = t.getDistanceAtXY( h1.x, h1.y, h1.z, duz, dv, h2.BzMid  );
	  } else {
	    err = t.getDistanceAtZ( h1.x, h1.y, h1.z, duz, dv, h2.BzMid  );
	  }
	  if( err!=0 ){
	    statCreateErrors++;
	    cout<<"Error by track chi2 calculation!!!: nhits="<<p.hits.size()<<" err= "<<err<<endl;
	    continue;
	  }	  
	  f[18] = duz;
	  f[19] = dv;
	  f[20] = t.getPtkG();
	  ntMCTracklets->Fill(f);

	} // ih3
      } // ih2
    } // ih1
  } // particles

  cout<<"--------\n"<<"tracklet cnstructor MC: "<<statCreateErrors<<" errors for "<<statCreateTries<<" created particles  ("<< 100.*statCreateErrors / statCreateTries <<" %)"<<endl; 

  statFile->Write();
}

