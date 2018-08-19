/*
 
  root -l statRecEff.root

  TNtuple *eff = (TNtuple*) gDirectory->FindObjectAny("eff");
  eff->Draw("r:z","w>0&&prim&&baseV==0&&pt<.3&&pt>=0.2&&hitRec!=1",""); 

  "part:nhits:pt:p:w" + ":prim:partR:partZ:nl:vol" + 
  ":layer:ihit:clust:x:y" + ":z:r:phi:fitpt:partRec" + 
  ":hitRec:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";

  // partRec: 0 or 1
  // hitRec: 
  // 0 not reconstructed
  // 1 reconstructed
  // 2 wrongly assigned
  // 3 assigned to ghost
  // 4 clone track

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


void TrackletConstructor::TrackletEfficiency( int useBaseLayers, const double ptCutMin, const double ptCutMax )
{ 
  /*
    0 - no base layers
    1 - 3 base layers
    2 - 2 base layers only tracks from the origin
  */

  if( !mTracker->mcFlag ) return;

  int iLayer1 = mCreationCuts.layerCuts[0].iLayer;
  int iLayer2 = mCreationCuts.layerCuts[1].iLayer;
  int iLayer3 = mCreationCuts.layerCuts[2].iLayer;
  
  if( useBaseLayers == 2  ){
    iLayer1 = -1;
    iLayer2 = mCreationPrimCuts.layerCuts[0].iLayer;
    iLayer3 = mCreationPrimCuts.layerCuts[1].iLayer;
  }

  
  static int nEvents = 0;
  static TFile *statFile = 0;
  static TNtuple *ntRecHits = 0;

  if( nEvents == 0 ){
    statFile = new TFile("statRecEff.root", "RECREATE");      
    string vars;
    vars = vars + 
      "part:nhits:pt:p:w" + ":prim:partR:partZ:nl:ihit" + 
      ":vol:layer:clust:x:y" + ":z:r:phi:fitpt:partRec" + 
      ":hitRec:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";

    ntRecHits = new TNtuple( "eff","reconstructed particle hits",vars.data() );
    // 0 not reconstructed
    // 1 reconstructed
    // 2 wrongly assigned
    // 3 assigned to ghost
    // 4 clone track
    ntRecHits->SetMarkerStyle(8);
    ntRecHits->SetMarkerSize(0.3);
  }
  nEvents++;


  vector<Tracklet> *tracks = &vTracks; 

  const int nTracks = tracks->size();
  const int nParticles = mTracker->mParticles.size();
  
  int partRecN[nParticles];  
  double partRecW[nParticles];  
  for( unsigned int i=0; i<nParticles; i++ ){
    partRecN[i] = 0;
    partRecW[i] = 0;
  }

  int recHits[mTracker->mHits.size()];
  for( int i=0; i<mTracker->mHits.size(); i++ ) recHits[i]=0;
     
  int trackRecPID[nTracks];
  int trackRecStatus[nTracks]; // 0 fake 1 rec 2 clone

  for( int i=0; i<nTracks; i++ ){
    trackRecPID[i] = -1;
    trackRecStatus[i] = 0;
  }


  static int nFakes[6]={0,0,0,0,0,0};
  static int nShort[6]={0,0,0,0,0,0};
  static double fakesW[6]={0,0,0,0,0,0};
  static double shortW[6]={0,0,0,0,0,0};

  for( int itr=0; itr<tracks->size(); itr++ ){
    Tracklet &tracklet = (*tracks)[itr];

    trackRecPID[itr] = -1;
    trackRecStatus[itr] = 0;

    // find the majority of mc indices

    int majorID = -1;
    int majorNhits = 0;
    double majorW=0;
    double totalW=0;
    {
      int trRecNHits[ nParticles ];
      double trRecW[ nParticles ];
      for( int i=0; i<nParticles; i++ ){
	trRecNHits[i]=0;
	trRecW[i] = 0;
      }

      for( int ih=0; ih<tracklet.vHits.size(); ih++ ){
	HitMC &mc = mTracker->mHitsMC[tracklet.vHits[ih]];
	if( mc.partID<0 ) majorNhits++;
	else {
	  trRecNHits[mc.partID]++;
	  trRecW[mc.partID]+= mc.w;
	}
	totalW+=mc.w;
      }
      
      for( int i=0; i<nParticles; i++ ){
	if( majorNhits < trRecNHits[i] ){
	  majorNhits = trRecNHits[i] ;
	  majorID = i;
	  majorW = trRecW[i];
	}
      }  
    }    
    
    int fakeType = 1;
    if( tracklet.nLayers>=5 ) fakeType = 5;
    else if( tracklet.nLayers==4 ) fakeType = 4;
    else if( tracklet.nLayers==3 ) fakeType = 3;
    else if( tracklet.nLayers==2 ) fakeType = 2;
    else  fakeType = 1;

    if( ( majorID<0 ) || 
	( majorNhits <= 0.5*tracklet.vHits.size() ) ){ // fake tracklet
      trackRecPID[itr] = -1; 
      trackRecStatus[itr] = 0;        
      fakesW[fakeType]+=totalW;
      nFakes[fakeType]++;
      fakesW[0]+=totalW;
      nFakes[0]++;
      continue;  
    }

    Particle &p = mTracker->mParticles[majorID];

    if( majorNhits <= 0.5 * p.hits.size() ){ // too short tracklet
       trackRecPID[itr] = majorID;
       trackRecStatus[itr] = 2;
       shortW[fakeType]+=totalW;      
       nShort[fakeType]++;
       shortW[0]+=totalW;      
       nShort[0]++;
      continue;
    }

    trackRecPID[itr] = majorID;
    trackRecStatus[itr] = 1;
 
    partRecN[majorID]++; 
    if( partRecW[majorID] < majorW ){
      partRecW[majorID] = majorW;    
    }
  } // tracklets

  

  // fill hit statuses for particle hits

  for( int itr=0; itr<tracks->size(); itr++ ){
    // 0 not reconstructed
    // 1 reconstructed
    // 2 wrongly assigned
    // 3 assigned to ghost
    // 4 clone track

    Tracklet &tracklet = (*tracks)[itr];
    
    int pid = trackRecPID[itr];
    int recStatus = trackRecStatus[itr]; // 0 fake 1 rec 2 clone
 
    for( int ih=0; ih<tracklet.vHits.size(); ih++ ){
      int hid = tracklet.vHits[ih];
      HitMC &mc = mTracker->mHitsMC[hid];
      
      if( recStatus == 0 ){// 0 fake track 
	recHits[hid] = 3;
	continue;
      }
      if( recStatus == 1 ){// 1 rec track
	if( mc.partID == pid ){
	  recHits[hid] = 1;
	} else {
	  recHits[hid] = 2;
	}
	continue;
      }
      if( recStatus == 2 ){// 2 clone track
	recHits[hid] = 4;
	continue;
      }
    }
  }


  // fill the ntuple

  for( unsigned int ipart=0; ipart<mTracker->mParticles.size(); ipart++ ){          
    Particle &p = mTracker->mParticles[ipart];         
    float f[100];
    for( int j=0; j<100; j++ ) f[j] = 0.;
    /*
      00 "part:nhits:pt:p:w" + ":prim:partR:partZ:nl:ihit" + 
      10 ":vol:layer:clust:x:y" + ":z:r:phi:fitpt:partRec" + 
      20 ":hitRec:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";
    */
    f[ 0] = ipart;
    f[ 1] = p.hits.size();
    f[ 2] = p.pt;
    f[ 3] = p.p;
    f[ 4] = p.w; 
    f[ 5] = p.prim;
    f[ 6] = p.r;
    f[ 7] = p.z;
    f[ 8] = p.nLayers;
    for( int ih=0; ih<p.hits.size(); ih++ ){
      int hid = p.hits[ih];
      Hit &h = mTracker->mHits[hid];
      f[ 9] = ih;
      f[10] = h.volume;
      f[11] = h.layer;
      f[12] = p.hitClusterIds[ih];
      f[13] = h.x;
      f[14] = h.y;
      f[15] = h.z;
      f[16] = h.r;
      f[17] = h.phi;
      f[18] = -1;
      f[19] = (partRecN[ipart]>0) ?1 :0;
      f[20] = recHits[hid];
      f[21] = p.baseV;
      f[22] = p.baseV0;
      f[23] = p.baseV1;
      f[24] = p.baseV2;
      f[25] = p.baseA0;
      f[26] = p.baseA1;
      f[27] = p.baseA2;
      f[28] = p.baseA3;
      
      ntRecHits->Fill(f);
    }
  }

  // fill the ntuple for fake hits
  {
    float f[100];
    for( int j=0; j<100; j++ ) f[j] = -1.;
    /*
      00 "part:nhits:pt:p:w" + ":prim:partR:partZ:nl:ihit" + 
      10 ":vol:layer:clust:x:y" + ":z:r:phi:fitpt:partRec" + 
      20 ":hitRec:bV:bV0:bV1:bV2" + ":bA0:bA1:bA2:bA3";
    */
    f[0] = -1;
    for( int hid=0; hid<mTracker->mHits.size(); hid++ ){
      Hit &h = mTracker->mHits[hid];
      f[10] = h.volume;
      f[11] = h.layer;
      f[13] = h.x;
      f[14] = h.y;
      f[15] = h.z;
      f[16] = h.r;
      f[17] = h.phi;
      f[20] = recHits[hid];
      ntRecHits->Fill(f);
    }
  }

  statFile->Write();
  

  // print efficiency
  

  static int particlesNCut = 0;
  static double particlesWCut = 0;
  static int particlesNRecCut = 0;
  static double particlesWRecCut = 0;
  static double particlesWRecFullCut = 0;
  
  for( unsigned int ipart=0; ipart<nParticles; ipart++ ){
    Particle &p = mTracker->mParticles[ipart];         
    //if( p.w<=0 ) continue;
    if( p.hits.size()<4 ) continue;
    //if( p.nLayers<3 ) continue;
    if(  p.pt<ptCutMin) continue;
    if(  ptCutMax > 0 && p.pt>=ptCutMax) continue;
    /*
      0 - no base layers
      1 - 3 base layers
      2 - 2 base layers only tracks from the origin
    */

    if( useBaseLayers==2){
      if( !p.baseV ) continue; //SG!!!
    } else {
      if( useBaseLayers==1 ){

	//if( p.baseV ) continue; //SG!!!
	//if( p.baseA0 ) continue; //SG
	int layersOK[3] = {0,0,0};
	for( int ih=0; ih<p.hits.size(); ih++ ){
	  int layer = mTracker->mHits[p.hits[ih]].layerID;
	  if( layer == iLayer1 ) layersOK[0] = 1;
	  if( layer == iLayer2 ) layersOK[1] = 1;
	  if( layer == iLayer3 ) layersOK[2] = 1;
	}
	int ok = ( layersOK[0] && layersOK[1] && layersOK[2] ) ;
	if( !ok ) continue;
      }
    }
    
    //if( !p.baseV ) continue; //SG!!!

    // cut on inner volumes
    //if( !p.prim ) continue; //SG
    //if( p.baseV>=0 ) continue; 
    //if( p.baseIC<0 ) continue; 
    //if( p.baseM>=0 ) continue; 
    //if( p.baseML>=0 ) continue; 
      //if( p.baseI<0 ) continue;      
    
    particlesNCut++;
    particlesWCut+=p.w;
    if( partRecN[ipart] > 0 ){
      particlesNRecCut++;
      particlesWRecCut+=partRecW[ipart];
      particlesWRecFullCut+=p.w;
    } else if(0) {
      cout<<"particle "<<ipart<<endl;
      p.Print();
      for( int ih=0; ih<p.hits.size(); ih++ ){
	Hit &h = mTracker->mHits[p.hits[ih]];
	h.Print();
      }
    }
  }

  static int statNTracks=0;
  
  statNTracks+=(*tracks).size();

  cout<<"------------- tracklet efficiency in "<<nEvents<<" events ---------"<<endl;
  cout<<"N tracklets "<<statNTracks<<endl;

  cout<<"Reconstructed particles: "<<100.*particlesNRecCut/particlesNCut      
      <<"% ("<<particlesNRecCut<<" out of "<<particlesNCut<<" particles )"<<endl;

  cout<<"Reconstructed w        : "<<100.*particlesWRecCut/particlesWCut      
      <<"% ("<<particlesWRecCut<<" out of "<<particlesWCut<<" weight )"<<endl;

  cout<<"Reconstructed hits w   : "<<100.*particlesWRecCut/particlesWRecFullCut
      <<"% ("<<particlesWRecCut<<" out of "<<particlesWRecFullCut<<" weight )"<<endl;
  
  cout<<"Fakes w                : "<<100.*fakesW[0]/particlesWCut<<"% (";
  for( int i=1; i<6; i++ ) cout <<i<<":"<<100.*fakesW[i]/particlesWCut<<"% ";
  cout<<") out of w "<<particlesWCut<<endl;

  cout<<"Shorts w               : "<<100.*shortW[0]/particlesWCut<<"% (";
  for( int i=1; i<6; i++ ) cout <<i<<":"<<100.*shortW[i]/particlesWCut<<"% ";
  cout<<") out of w "<<particlesWCut<<endl;


  cout<<"Fakes n                : "<<nFakes[0]<<" (";
  for( int i=1; i<6; i++ ) cout <<i<<":"<<nFakes[i]<<" ";
  cout<<") out of "<<statNTracks<<" tracks:"<<endl;

  cout<<"Shorts n               : "<<nShort[0]<<" (";
  for( int i=1; i<6; i++ ) cout <<i<<":"<<nShort[i]<<" ";
  cout<<") out of "<<statNTracks<<" tracks:"<<endl;


}


