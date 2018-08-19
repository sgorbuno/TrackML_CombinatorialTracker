#include "Tracker.h"
#include "util.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include <iostream>
#include "TrackModelPhysical.h"
#include "TrackletConstructor.h"
#include "TMath.h"


int layerNHits[Geo::NLayers];

bool checkLayers( const int *baseLayers, int nBaseLayers )
{
  int icross=0;
  for( int il=0; il<nBaseLayers; il++ ){
    if( layerNHits[baseLayers[il]]>0 ) icross++;    
  }
  return ( icross==nBaseLayers );
}


void Tracker::readEvent( const char *directory, int event, bool loadMC )
{
  TString filePrefix;
  filePrefix.Form("%sevent%09d",directory,event);

  mcFlag = loadMC;

  mHits.clear();
  //mHitRecInfo.clear();
  //mTracks.clear();
  mHitsMC.clear();
  mParticles.clear();
  mHits.reserve(150000); 

  // ===== load hits
  { 
    TString fname = filePrefix+"-hits.csv";
    ifstream in(fname.Data());    
    if( !in.is_open() ){
      cout<<"Event "<<event<<" does not exist!!"<<endl;
      exit(0);
    }
    char tmpLine[256];
    in.getline(tmpLine,256);
    cout<<tmpLine<<endl;
    while (1) {    
      double h[7]; //file line: id:x:y:z:volume:layer:module
      if( !readLine(in,h,7) ) break;       
      if( h[0]-1 != mHits.size() ){
	cout<<"Hit index is wrong: "<<h[0]<<endl;
	exit(0);
      }      
      Hit hit;
      hit.x = h[1]*0.1; // convert to [cm]
      hit.y = h[2]*0.1; // convert to [cm]
      hit.z = h[3]*0.1; // convert to [cm]      
      hit.r = sqrt( hit.x*hit.x + hit.y*hit.y );
      hit.phi = atan2( hit.y, hit.x );
      hit.module = h[6];
      hit.layer = ( (int) h[5])/2 -1;
      switch( (int) h[4] )
	{
	case  8: hit.volume = 0; break;
	case  7: hit.volume = 1; hit.layer = 6-hit.layer; break;
	case  9: hit.volume = 2; break;
	case 13: hit.volume = 3; break;
	case 12: hit.volume = 4; hit.layer = 5-hit.layer; break;
	case 14: hit.volume = 5; break;
	case 17: hit.volume = 6; break;
	case 16: hit.volume = 7; hit.layer = 5-hit.layer; break;
	case 18: hit.volume = 8; break;
	default:	  
	  cout<<"Unknown detector volume: "<< (int) h[4] << endl;
	  exit(0);	
	};     

      if( hit.layer<0 || hit.layer>= Geo::volumes[hit.volume].nLayers ){
	cout<<"Unknown detector layer: "<<hit.layer<<endl;
	exit(0);	
      }
      hit.layerID = Geo::volumes[hit.volume].layerIDs[hit.layer];
      Layer &layer = Geo::layers[hit.layerID];
      if( Geo::volumes[hit.volume].type==0 ){
	hit.t = hit.z / hit.r * layer.r;
      } else {
	hit.t = hit.r / hit.z * layer.z ;
      }
      hit.isUsed=0;
      hit.trackID=0;
      hit.BzFwd = Geo::OriginBzkG;
      hit.BzMid = Geo::OriginBzkG;
      hit.BzBck = Geo::OriginBzkG;
      if( 1|| ( hit.volume==4 || hit.volume==5 || hit.volume==7 || hit.volume==8 ) ){
	hit.BzFwd = Geo::layers[hit.layerID].getFieldFwd(hit.phi, hit.t);
 	hit.BzMid = Geo::layers[hit.layerID].getFieldMid(hit.phi, hit.t);
	hit.BzBck = Geo::layers[hit.layerID].getFieldBck(hit.phi, hit.t);
       }
      mHits.push_back(hit);
    }
    cout<<" loaded "<<mHits.size()<<" hits "<<endl;
    in.close();    

  } // load hits

 
  if( !loadMC ){
    for( int ilayer=0; ilayer<Geo::NLayers; ilayer++ ){
      FitLayer &layer = mFitLayers[ilayer];
      layer.mTracker = this;
      layer.Create( ilayer, 1, 5., 5. );
    }
    return;
  }
  
  
  // create particle ID->index map
  std::map<long unsigned int,int> partIDmap;   

  // ========= load particles with reindexing
  {
    mParticles.clear();          

    TString fname = filePrefix+"-particles.csv";
    ifstream in(fname.Data());
    if( !in.is_open() ){
      cout<<"Particle file for event "<<event<<" does not exist!!"<<endl;
      exit(0);
    }
    char tmpLine[256];
    in.getline(tmpLine,256);
    cout<<tmpLine<<endl;
    while (1) {    
      double f[9]; //  particle_id,vx,vy,vz,px,py,pz,q,nhits
      if( !readLine(in,f,9) ) break;
      int nhits = (int) f[8];
      if( nhits==0 ) continue; // no hits belong to the particle		
      Particle p( nhits );
      p.x = f[1]*0.1; // [cm]
      p.y = f[2]*0.1; // [cm]
      p.z = f[3]*0.1; // [cm]
      p.r = sqrt(p.x*p.x+p.y*p.y);
      p.px = f[4];
      p.py = f[5];
      p.pz = f[6];
      p.q  = f[7];
      p.xl = 0;
      p.yl = 0;
      p.zl = 0;
      p.rl = 0;
      p.pt = sqrt(f[4]*f[4] + f[5]*f[5]);
      p.p  = sqrt(f[4]*f[4] + f[5]*f[5] + f[6]*f[6]);
      p.prim = fabs(p.z)<1.2 && p.r<0.05;
      p.w = 0;      
      partIDmap[ (long unsigned int) f[0] ] = mParticles.size();
      mParticles.push_back(p);
    }
    cout <<" loaded "<<mParticles.size() <<"particles in event "<<event<<endl;    
    in.close();
  } // particles


  { // ============  read  mc truth

    mHitsMC.clear();
    mHitsMC.reserve(150000);

    TString fname = filePrefix+"-truth.csv";
    ifstream in(fname.Data());    
    if( !in.is_open() ){
      cout<<"Truth file for event "<<event<<" does not exist!!"<<endl;
      exit(0);
    }

    char tmpLine[256];    
    in.getline(tmpLine,256);
    cout<<tmpLine<<endl;

    while (1) {    
      double mc[9]; //mc: hit_id,particle_id,tx,ty,tz,tpx,tpy,tpz,weight     
      if( !readLine(in,mc,9) ) break;
      if( mc[0]-1 != mHitsMC.size() ){
	cout<<"MC hit index is wrong: "<<mc[0]<<endl;
	exit(0);
      } 
      HitMC hitmc;
      Hit &hit = mHits[mHitsMC.size()];
 
      hitmc.hitID = mHitsMC.size();
      hitmc.x = mc[2]*0.1; // convert to [cm]
      hitmc.y = mc[3]*0.1; // convert to [cm]
      hitmc.z = mc[4]*0.1; // convert to [cm]      
      hitmc.partID = -1;
      hitmc.w = mc[8]; // weight 
      hitmc.px = mc[5];
      hitmc.py = mc[6];
      hitmc.pz = mc[7];
      hitmc.pt = sqrt(hitmc.px*hitmc.px + hitmc.py*hitmc.py); // pt
      hitmc.dOrigin2 = -1.; // squared distance to particle origin 
      hitmc.p = sqrt( hitmc.px*hitmc.px + hitmc.py*hitmc.py + hitmc.pz*hitmc.pz ); // momentum  for sorting
      hitmc.q = -1000;
      hitmc.r = sqrt( hitmc.x*hitmc.x + hitmc.y*hitmc.y );
      hitmc.phi = atan2( hitmc.y, hitmc.x );     
      if( Geo::volumes[hit.volume].type==0 ) hitmc.t = hitmc.z;
      else hitmc.t = hitmc.r;


      // find mapped particle id
      long unsigned int id = mc[1];
      int newID = 0;
      if( id==0 ){ // hit is not associated to any particle 
	mHitsMC.push_back(hitmc);
	continue;     	
      }
      std::map<long unsigned int, int>::iterator it = partIDmap.find(id);
      if( it==partIDmap.end() || it->first!=id){
	cout<<"Mapped particle ID is wrong!!!"<<endl;
	cout<<"ID= "<<id<<" hit "<<id<<" iterator at ID "<<it->first<<endl;
	exit(0);
      }
      newID = it->second;
      if( newID < 0 || newID>= mParticles.size() ){
	cout<<"Mapped particle ID is wrong!!!"<<endl;
	cout<<"ID= "<<id<<" new ID "<<newID<<endl;
	exit(0);
      }
      hitmc.partID = newID;     
      Particle &p = mParticles[newID];

      double dx = hitmc.x - p.x;
      double dy = hitmc.y - p.y;
      double dz = hitmc.z - p.z;
      hitmc.dOrigin2 = dx*dx + dy*dy + dz*dz;

      p.w += hitmc.w;
      hitmc.q = p.q;

      p.hits.push_back( hitmc.hitID );
      mHitsMC.push_back(hitmc);           
    }
    cout << " read "<<mHitsMC.size() << " mc hits for event "<<event<<endl;
    in.close();      
    if( mHitsMC.size() != mHits.size() ){
      cout<<"number of MC hits is not equal to the number of hits"<<endl;
      exit(0);
    }
  } // read mc hits
    
  
  // sort particle hits

  for( unsigned int ipart=0; ipart<mParticles.size(); ipart++ ){      
    Particle &p = mParticles[ipart];
    if( p.hits.size()<1 ) continue;
    std::vector<HitMC> v;
    for( int i=0; i<p.hits.size(); i++ ){
      v.push_back( mHitsMC[p.hits[i]] );
    }
    std::sort(v.begin(), v.end() );
    for( int i=0; i<p.hits.size(); i++ ){
      p.hits[i] = v[i].hitID;
    }
    p.hitClusterIds.resize(p.hits.size());
    int il=-1;
    int lastLid=-1;
    for( int i=0; i<p.hits.size(); i++ ){
      if( mHits[p.hits[i]].layerID != lastLid ){
	il++;
	lastLid = mHits[p.hits[i]].layerID;
      }
      p.hitClusterIds[i] = il;
    } 

    {
      int ihlast = p.hits[p.hits.size()-1];
      HitMC &mc = mHitsMC[ihlast];      
      p.xl = mc.x;
      p.yl = mc.y;
      p.zl = mc.z;
      p.rl = sqrt(mc.x*mc.x+mc.y*mc.y);
    }
    if( 0 && ipart==467 ){
      for( int i=0; i<p.hits.size(); i++ ){
	Hit &h = mHits[p.hits[i]];      
	h.Print();
      }
      for( int i=0; i<p.hits.size(); i++ ){
	HitMC &mc = mHitsMC[p.hits[i]];      
	mc.Print();	
      }
    }
  }
  //exit(0);
  

  for( unsigned int ipart=0; ipart<mParticles.size(); ipart++ ){      
    Particle &p = mParticles[ipart];
    for( int i=0; i<Geo::NLayers; i++ ) layerNHits[i]=0;
    for( int i=0; i<p.hits.size(); i++ ){
      layerNHits[mHits[p.hits[i]].layerID]++;      
    }
    p.nLayers = 0;
    for( int i=0; i<Geo::NLayers; i++ ){
      if( layerNHits[i]>0 ) p.nLayers ++;
    }

    p.baseV0 = 0;
    p.baseV1 = 0;
    p.baseV2 = 0;
    p.baseV  = 0;
    if( p.prim ){
      int baseV0[2] ={ Geo::getLayerID(0,0), Geo::getLayerID(0,1) };
      int baseV1[2] ={ Geo::getLayerID(1,0), Geo::getLayerID(1,1) };
      int baseV2[2] ={ Geo::getLayerID(2,0), Geo::getLayerID(2,1) };
      p.baseV0 = checkLayers(baseV0,2);
      p.baseV1 = checkLayers(baseV1,2);
      p.baseV2 = checkLayers(baseV2,2);
      p.baseV  =  ( p.baseV0 || p.baseV1 || p.baseV2 );
    }

    int baseA0[3] ={ Geo::getLayerID(0,0), Geo::getLayerID(0,1), Geo::getLayerID(0,2) };
    p.baseA0 = checkLayers(baseA0,3);

    int baseA1[3] ={ Geo::getLayerID(1,0), Geo::getLayerID(1,1), Geo::getLayerID(1,2) };
    p.baseA1 = checkLayers(baseA1,3);

    int baseA2[3] ={ Geo::getLayerID(2,2), Geo::getLayerID(2,3), Geo::getLayerID(2,4) };
    p.baseA2 = checkLayers(baseA2,3);

    int baseA3[3] ={ Geo::getLayerID(3,0), Geo::getLayerID(3,1), Geo::getLayerID(3,2) };
    p.baseA3 = checkLayers(baseA3,3);

  } // particles


 
  // create fit layers
  
  for( int ilayer=0; ilayer<Geo::NLayers; ilayer++ ){
    FitLayer &layer = mFitLayers[ilayer];
    layer.mTracker = this;
    layer.Create( ilayer, 1, 5., 5. );
  }

}
