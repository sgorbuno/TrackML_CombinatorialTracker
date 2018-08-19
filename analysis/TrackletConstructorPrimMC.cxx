/*
  root -l statMCTracklets.root
  TNtuple *seeds = (TNtuple*) gDirectory->FindObjectAny("mcTracklets");  
  
  "part:nhits:pt:p:w" + ":nl:prim:layer:dupl:x" + 
  ":y:z:r:phi:t" + ":dPL:dTL:dTH:duz:dv" +
  ":fitpt:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";

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


void TrackletConstructor::ReconstructPrimMC( const CreationPrimCuts &crCuts, const ProlongationCuts &prCuts )
{  
  mCreationPrimCuts = crCuts;
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

    ntMCTracklets = new TNtuple ("mcTracklets","mc tracklet layers", vars.data() );
    ntMCTracklets->SetMarkerStyle(8);
    ntMCTracklets->SetMarkerSize(0.3);
  }
  nEvents++;

  static int statCreateTries=0;
  static int statCreateErrors=0;

  int ilayer1 = crCuts.layerCuts[0].iLayer;
  int ilayer2 = crCuts.layerCuts[1].iLayer;
  FitLayer layer1;
  layer1.mTracker = mTracker;
  layer1.Create( ilayer1, 0, 0.2, 4. );
  FitLayer layer2;
  layer2.mTracker = mTracker;
  layer2.Create( ilayer2, 0, 0.2, 4. );

  
  for( unsigned int ipart=0; ipart<mTracker->mParticles.size(); ipart++ ){      

    Particle &p = mTracker->mParticles[ipart];     

    // check the number of crossed layers
    int partNlayers = 0;
    {
      int present[2] = {0,0};
      for( int ih = 0; ih<p.hits.size(); ih++ ){
	Hit &h = mTracker->mHits[p.hits[ih]];
	if( h.layerID == ilayer1) present[0] = 1;
	if( h.layerID == ilayer2) present[1] = 1;
      }       
      if( !present[0] ) partNlayers = 0;
      else if( !present[1] ) partNlayers = 1;
      else partNlayers = 2;
    }
    
    if( partNlayers < 2 ) continue;
    /*
      00 "part:nhits:pt:p:w" + ":nl:prim:layer:dupl:x" + 
      10 ":y:z:r:phi:t" + ":dPL:dTL:dTH:duz:dv" +
      20 ":fitpt:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";
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

    int duplicate[2] = {-1,-1};

    duplicate[0] = -1;
    
    for( int ih1 = 0; ih1<p.hits.size(); ih1++ ){
      
      Hit h1 = mTracker->mHits[p.hits[ih1]];
      if( h1.layerID != ilayer1 ) continue;      

      duplicate[0]++;       
      /*
	00 "part:nhits:pt:p:w" + ":nl:prim:layer:dupl:x" + 
	10 ":y:z:r:phi:t" + ":dPL:dTL:dTH:duz:dv" +
	20 ":fitpt:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";
      */
      f[ 7] = 1;
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

      double extr2Tslope = 0;
      if( layer2.mType==0 ){	  
	extr2Tslope =  h1.z/h1.r;
      } else {
	extr2Tslope =  h1.r/h1.z;
      }

      duplicate[1] = -1;
      for( int ih2 = 0; ih2<p.hits.size(); ih2++ ){
	Hit &h2 = mTracker->mHits[p.hits[ih2]];
	if( h2.layerID != ilayer2 ) continue;			
	duplicate[1]++;
         /*
	   00 "part:nhits:pt:p:w" + ":nl:prim:layer:dupl:x" + 
	   10 ":y:z:r:phi:t" + ":dPL:dTL:dTH:duz:dv" +
	   20 ":fitpt:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";
	 */
	f[ 7] = 2;
	f[ 8] = duplicate[1];
	f[ 9] = h2.x;
	f[10] = h2.y;
	f[11] = h2.z;
	f[12] = h2.r;
	f[13] = h2.phi;
	f[14] = h2.t;
	f[15] = h2.phi - h1.phi; // extrapolated to layer phi
	double extrTL = 0, extrTH = 0;	 
	if( layer2.mType==0 ){
	  extrTL = extr2Tslope*(layer2.mR );
	  extrTH = extr2Tslope*(  h2.r  );
	  f[16] = h2.t - extrTL;// extrapolated to layer T	
	  f[17] = h2.z - extrTH;// extrapolated to hit T
	} else {
	  extrTL = extr2Tslope*(layer2.mZ );
	  extrTH = extr2Tslope*(    h2.z  );
	  f[16] = h2.t - extrTL;// extrapolated to layer T	
	  f[17] = h2.r - extrTH;// extrapolated to hit T
	}
	
	// create tracklet	  
	statCreateTries++;
	TrackModelPhysical t;
	int err = t.createXYLast( 0., 0., 0.,
				  h1.x, h1.y, h1.z,			      
				  h2.x, h2.y, h2.z , h1.BzMid );
	
	if( err!=0 ){
	  statCreateErrors++;
	  cout<<"Error by tracklet creation!!!: prim "<<p.prim<<" nhits="<<p.hits.size()<<" err= "<<err<<endl;
	  cout<<h1.x<<" "<<h1.y<<" "<<h1.z<<" "<<endl;
	  cout<<h2.x<<" "<<h2.y<<" "<<h2.z<<" "<<endl;	 
	  continue;
	}
	
	double dv=0, duz=0;	  
	err = t.getDistanceAtXY( 0, 0, 0, duz, dv, h1.BzMid  );
	if( err!=0 ){
	  statCreateErrors++;
	  cout<<"Error by track chi2 calculation!!!: nhits="<<p.hits.size()<<" err= "<<err<<endl;
	  continue;
	}	  
	f[18] = duz;
	f[19] = dv;
	f[20] = t.getPtkG();

	ntMCTracklets->Fill(f);
      } // ih2        
    } // ih1
  } // particles

  cout<<"--------\n"<<"tracklet cnstructor MC: "<<statCreateErrors<<" errors for "<<statCreateTries<<" created particles  ("<< 100.*statCreateErrors / statCreateTries <<" %)"<<endl; 

  statFile->Write();
}

