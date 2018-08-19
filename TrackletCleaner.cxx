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


struct SortHit
{
  int isUsed;
  int layer; 
};

struct SortTracklet
{
  bool isAlive;
  double chi2; 
  int nActiveHits;
  int nActiveLayers;
  int ptRegion;
};


void TrackletConstructor::CleanTracklets( int minTrackNLayers )
{  
  cout<<"Sorting "<<vTracklets.size()<<" tracklets... "<<endl;

  constexpr int nLayers = Geo::NLayers;

  if( vTracklets.size()==0 ) return;

  vector<Tracklet> vTr = move(vTracklets);
  int nHits = mTracker->mHits.size();
  int nTracklets = vTr.size();

  SortHit *sHits = new SortHit[nHits];
  SortTracklet *sTracklets = new SortTracklet[nTracklets];
 
  for( int i=0; i<nHits; i++ ){
    Hit &h = mTracker->mHits[i];
    SortHit &sh = sHits[i];
    sh.isUsed = h.isUsed;
    sh.layer = h.layerID;
  }

  int nTrHitsPerLayer[nLayers];

  for( int itr=0; itr<nTracklets; itr++ ){
     Tracklet &t = vTr[itr];
     SortTracklet &st = sTracklets[itr];
     st.isAlive = 1;    
     st.chi2 = t.chi2/t.vHits.size();
     double ptstep = 1./0.05;
     st.ptRegion = (int) t.pt*ptstep;
     if( st.ptRegion > (int) 1.*ptstep ) st.ptRegion = (int) 1.*ptstep;

     // number of layers 
     for( int il=0; il<nLayers; il++ ){
       nTrHitsPerLayer[il]=0;
     }
     for( int ih=0; ih<t.vHits.size(); ih++ ){
       SortHit &sh = sHits[t.vHits[ih]];
       nTrHitsPerLayer[sh.layer]++;
     }
     t.nLayers=0;
     for( int il=0; il<nLayers; il++ ){
       if( nTrHitsPerLayer[il]>0 ) t.nLayers++;
     }   
     st.nActiveLayers = t.nLayers;
     st.nActiveHits = t.vHits.size();
  }

  int firstAlive=0;
  int lastAlive=nTracklets-1;
  do{
    
    while( firstAlive<=lastAlive && !sTracklets[firstAlive].isAlive ) firstAlive++;
    while( firstAlive<=lastAlive && !sTracklets[lastAlive].isAlive ) lastAlive--;

    int bestTrId = -1;
    int bestTrNHits = 0;
    int bestTrPtRegion=0;
    int bestTrVtx = 0;
    double bestTrChi2 = 1.e10;

    for( int itr=firstAlive; itr<=lastAlive; itr++ ){
      Tracklet &t = vTr[itr];
      SortTracklet &st = sTracklets[itr];
      if( !st.isAlive ) continue;      
      
      if( st.nActiveHits < bestTrNHits ) continue;
      
      st.nActiveHits = 0;
      for( int ih=0; ih<t.vHits.size(); ih++ ){
	SortHit &sh = sHits[t.vHits[ih]];
	if( sh.isUsed ) continue;
	st.nActiveHits++;	
      }
      if( st.nActiveHits < 2 ){
	st.isAlive=0;
	continue;
      }
      if( bestTrNHits > st.nActiveHits ) continue;
      if( bestTrNHits==st.nActiveHits ){
	//if( bestTrPtRegion > st.ptRegion  ) continue; 
	//if( bestTrPtRegion >= st.ptRegion  && bestTrChi2 < st.chi2) continue; 
	if( bestTrVtx && !t.vtxTrack ) continue;
	if(  bestTrChi2 < st.chi2) continue; 
      }             
      
      // check the number of missed layers 
      for( int il=0; il<nLayers; il++ ){
	nTrHitsPerLayer[il]=0;
      }
      for( int ih=0; ih<t.vHits.size(); ih++ ){
	SortHit &sh = sHits[t.vHits[ih]];
	if( sh.isUsed) continue;
	nTrHitsPerLayer[sh.layer]++;
      }
      st.nActiveLayers = 0;
      for( int il=0; il<nLayers; il++ ){
	if( nTrHitsPerLayer[il]>0 ) st.nActiveLayers ++;
      }   
      
      if( st.nActiveLayers < t.nLayers - 3 ){
	st.isAlive=0;
	continue;
      }
      
      bestTrId = itr;
      bestTrNHits = st.nActiveHits;
      bestTrChi2 = st.chi2;
      bestTrPtRegion = st.ptRegion;
      bestTrVtx = t.vtxTrack;
    }
    
    if( bestTrNHits < 2 ) break;

    {
      Tracklet &t = vTr[bestTrId];
      SortTracklet &st = sTracklets[bestTrId];
      Tracklet tt = t;
      tt.vHits.clear();
      tt.nLayers = st.nActiveLayers;
      st.isAlive = 0;
      bool writeToTrack = ( tt.nLayers >= minTrackNLayers );
      for( int ih=0; ih<t.vHits.size(); ih++ ){
	SortHit &sh = sHits[t.vHits[ih]];
	Hit &h = mTracker->mHits[t.vHits[ih]];
	if( h.isUsed && !sh.isUsed ){
	  cout<<"sorting: wrong indexing"<<endl;
	  exit(1);
	}
	if( sh.isUsed ) continue;
	tt.vHits.push_back(t.vHits[ih]);
	sh.isUsed = 1;
	if( writeToTrack ){
	  if( h.isUsed!=0 ){
	    cout<<"reconstructed already used hit!!! hit n "<<t.vHits[ih]<<endl;
	    h.Print();
	    exit(1);
	  }
	  h.isUsed = 1;
	}
      }

      // write to track array or back to the tracklets array
      if( writeToTrack ){ 
	vTracks.push_back(tt);
      } else {
	vTracklets.push_back(tt);
      }
    }
  } while( 1);


  delete[] sHits;
  delete[] sTracklets;

  cout<<"...sorting"<<endl;

}
