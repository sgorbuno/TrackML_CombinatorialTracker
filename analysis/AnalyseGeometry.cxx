/*
  root -l statGeometry.root

  TNtuple *geo = (TNtuple*) gDirectory->FindObjectAny("geo");  
  
  // 0:part :nhits :pt :p :q 5:w :nl :prim :vol :lay 
  //10:layID :mx :my :mz :mr 15:mphi :mt :mpt :mp  :hx 
  //20:hy :hz :hr :hphi :ht 25:Bz0 :Bz1 :Bz2 :fitBz0 :fitBz1
  //30:fitBz2

*/

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
#include "Geo.h"
#include "Cuts.h"
#include <iostream>
#include "PolynomFit.h"
#include "TMath.h"

using namespace std;

  constexpr int NFieldPar=9;

class StatLayer
{
public:
  int id=0;
  int vol=0;
  int layer=0;
  double tMin=1.e10;
  double tMax=1.e-10;  
  double r=0;
  double z=0;
  long int nEntriesRZ=0;

  PolynomFit fit[3]; // fwd, mid, bck

  
  double cField[3][NFieldPar]; // Bz field in kGaus (t,phi) for fwd, fit, bckwd
  double fieldSigma[3]={0,0,0};

  double getField( int region, double phi, double t){
    double t2 = t*t;
    double *c=cField[region];
    return 
      (c[0] + c[1]*t + c[2]*t2) + 
      (c[3] + c[4]*t + c[5]*t2)*sin(phi) + 
      (c[6] + c[7]*t + c[8]*t2)*cos(phi);
  }

  double fieldS2[3]={0,0,0};
  long int fieldS2Entries[3]={0,0,0};
};


StatLayer statLayers[Geo::NLayers];

void resetLayerFit()
{
 for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    for( int ir=0; ir<3; ir++ ){
      sl.fit[ir].Reset( NFieldPar );
    }
 }
}

void initLayers()
{
  for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    Layer &l = Geo::layers[il];
    sl.id = il;
    sl.vol = l.volume;
    sl.layer = l.layer;
    for( int ir=0; ir<3; ir++ ){
      for( int i=0; i<NFieldPar; i++){
	sl.cField[ir][i]=0;
      }    
      sl.fieldSigma[ir]=0;
      sl.fieldS2[ir]=0;
      sl.fieldS2Entries[ir]=0;
   }
  }
  resetLayerFit();
}

void writeLayerSizes()
{
  ofstream out("geoLayerSizes.txt");    
  if( !out.is_open() ){
    cout<<"analyseGeometry:: Can not open output file"<<endl;
    exit(0);
  }
  for( int i=0; i<Geo::NLayers; i++){
    StatLayer &l = statLayers[i];
    if( l.nEntriesRZ>0 ){
      l.r/=l.nEntriesRZ;
      l.z/=l.nEntriesRZ;
      l.nEntriesRZ=1;
    }
    out<<i<<" "<<l.r<<" "<<l.z<<" "<<l.tMin<<" "<<l.tMax<<endl;
  }
  out.close();
}

void readLayerSizes()
{
  ifstream in("geoLayerSizes.txt");    
  if( !in.is_open() ){
    cout<<"analyseGeometry:: Can not open input file"<<endl;
    exit(0);
  }
  for( int i=0; i<Geo::NLayers; i++){
    StatLayer &l = statLayers[i];
    int j;
    in>>j>>l.r>>l.z>>l.tMin>>l.tMax;
    if( j!=i ){
      cout<<"geo file broken"<<endl;
      exit(1);
    }
    l.nEntriesRZ = 1;
  }
  in.close();
}

void updateLayerSizesRZ( Tracker *mTracker)
{
  for( int ih=0; ih<mTracker->mHits.size(); ih++ ){
    const Hit &h = mTracker->mHits[ih];    
    StatLayer &sl = statLayers[h.layerID];
    sl.r+= h.r;
    sl.z+= h.z;
    sl.nEntriesRZ++;
  }
}

void updateLayerSizesT( Tracker *mTracker)
{
  for( int ih=0; ih<mTracker->mHits.size(); ih++ ){
    const Hit &h = mTracker->mHits[ih];    
    Layer &l=Geo::layers[h.layerID];
    StatLayer &sl = statLayers[h.layerID];

    double t = 0;
    if( l.type==0 ) t = h.z;
    else  t = h.r;
    if( sl.tMin > t ) sl.tMin = t;
    if( sl.tMax < t ) sl.tMax = t;    

    t = 0;
    if( l.type==0 ) t = h.z/h.r*sl.r;
    else  t = h.r/h.z*sl.z;
    if( sl.tMin > t ) sl.tMin = t;
    if( sl.tMax < t ) sl.tMax = t;    
    
  }
}


void writeLayerField()
{
  ofstream out("geoLayerField.txt");    
  if( !out.is_open() ){
    cout<<"analyseGeometry:: Can not open output file"<<endl;
    exit(0);
  }

  for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    for( int ir=0; ir<3; ir++){
      sl.fieldSigma[ir] = -1;
      if( sl.fieldS2Entries[ir]>0 ) sl.fieldSigma[ir] = sqrt(sl.fieldS2[ir]/sl.fieldS2Entries[ir]);
      out<<il;
      for( int i=0; i<NFieldPar; i++)  out<<" "<<sl.cField[ir][i];
      out<<" "<<sl.fieldSigma[ir]<<endl;
    }
  }
  out.close();
}

void readLayerField()
{
  ifstream in("geoLayerField.txt");    
  if( !in.is_open() ){
    cout<<"analyseGeometry:: Can not open input file"<<endl;
    exit(0);
  }
  for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    for( int ir=0; ir<3; ir++){
      int jl;
      in>>jl;
      if( jl!=il ){
	cout<<"geo field file broken"<<endl;
	exit(1);
      }
      for( int i=0; i<NFieldPar; i++ ) in>>sl.cField[ir][i];
      in>>sl.fieldSigma[ir];
      sl.fieldS2[ir]=0;
      sl.fieldS2Entries[ir]=0;
    }
  }
  in.close();
}


void fitLayerField()
{
  for( int il=0; il<Geo::NLayers; il++){
    StatLayer &sl = statLayers[il];
    for( int ir=0; ir<3; ir++ ){
      int err = sl.fit[ir].Fit(sl.cField[ir]);
      if( err!=0 ){
	cout<<"Can not fit the field!!!! with "<<sl.fieldS2Entries[ir]<<" measurements, error "<<err<<endl;
	sl.fieldSigma[ir] = -1;
      }
    }
  }
}

void updateLayerField( Tracker *mTracker, vector<double> hitField[3], double cut )
{
  for( int ih=0; ih<mTracker->mHits.size(); ih++ ){
    const Hit &h = mTracker->mHits[ih];
    for( int ir=0; ir<3; ir++){
      double mesBz = hitField[ir][ih];
      if( mesBz<-100. ) continue;
      //Layer &l=Geo::layers[h.layerID];
      StatLayer &sl = statLayers[h.layerID];
      double polBz = sl.getField( ir, h.phi, h.t );
      double d = mesBz - polBz;
      bool specialRegion = (h.layerID==35) || ((h.layerID==34)&&ir!=2);
      if( !specialRegion && cut>0. && fabs(d)>cut*sl.fieldSigma[ir] ) continue;
      sl.fieldS2[ir]+=d*d;
      sl.fieldS2Entries[ir]++;
      double s = sin(h.phi);
      double c = cos(h.phi);
      double t = h.t;
      double t2 = t*t;
      double f[NFieldPar] = {1, t, t2, s, t*s, t2*s, c, t*c, t2*c};
      sl.fit[ir].AddMeasurement( f, mesBz );    
    }
  }
}



void Tracker::analyzeGeometry( bool endOfData )
{ 
  cout<<"\n-----------------\n analyse geometry  .."<<endl;

  static int nEvents = 0;
  static TFile *statFile = 0;
  static TNtuple *ntGeo = 0;
  
  if( nEvents == 0 ){
    statFile = new TFile("statGeometry.root", "RECREATE");  

    ntGeo = new TNtuple
      ("geo","geometry",
       "part:nhits:pt:p:q:w:nl:prim:vol:lay:layID:mx:my:mz:mr:mphi:mt:mpt:mp:hx:hy:hz:hr:hphi:ht:Bz0:Bz1:Bz2:fitBz0:fitBz1:fitBz2");
    ntGeo->SetMarkerStyle(8);    
    ntGeo->SetMarkerSize(0.3);

    initLayers();
    readLayerSizes();
    readLayerField();
  }
  nEvents++;

  if( endOfData ){
    statFile->Write();
    statFile->Close();
    nEvents=0;
    //writeLayerSizes();
    fitLayerField();
    //writeLayerField();
    return;
  }

  //updateLayerSizesT(this);  
  //return;

  // find magnetic field value at each hit

  int nHits = mHits.size();
  vector<double> vBz[3];
  for( int ir=0; ir<3; ir++){
    vBz[ir].resize(nHits);
    for( int ih=0; ih<nHits; ih++) vBz[ir][ih] =-1000;
  }

  for( unsigned int ipart=0; ipart<mParticles.size(); ipart++ ){          
    Particle &p = mParticles[ipart];             
    if( !p.prim ) continue;

    for( int ih=0; ih<p.hits.size(); ih++ ){

      int hid = p.hits[ih];
      Hit &h = mHits[hid];
	   
      // search for the next hit
      int iFwd = -1;
      for( int ih1=ih+1; ih1<p.hits.size(); ih1++ ){
	int hid1 = p.hits[ih1];
	Hit &h1 = mHits[hid1];
	if( (h1.layerID > h.layerID) && (fabs(h1.z-h.z)<200.) ){	  
	  iFwd = ih1;
	  break;
	}
      }
      int iBck=-1;
      { // search backwards
	for( int ih1=ih-1; ih1>=0; ih1-- ){
	  int hid1 = p.hits[ih1];
	  Hit &h1 = mHits[hid1];
	  if( (h1.layerID < h.layerID) && (fabs(h1.z-h.z)<200.)){
	    iBck = ih1;
	    break;
	  }
	}
      }
      
      HitMC mc = mHitsMC[hid];

      if( mc.pt>1. || mc.pt<.2 ) continue;

      HitMC mcFwd;
      if( iFwd>=0 ) mcFwd = mHitsMC[p.hits[iFwd]];
      HitMC mcBck;
      mcBck.x=0;
      mcBck.y=0;
      mcBck.z=0;      
      if( iBck>=0 ) mcBck = mHitsMC[p.hits[iBck]];
      
      if( mc.partID!=ipart || mc.hitID!=hid || mc.q!=p.q ||
	  (iFwd>=0&&mcFwd.partID!=ipart) || (iBck>=0&&mcBck.partID!=ipart)  ){
	cout<<"geometry: hit indexation is broken!!"<<endl;
	exit(1);
      }
      
      // fit forward field
      if( iFwd>=0 ){
	if( mc.px*mcFwd.px+mc.py*mcFwd.py < 0 ){
	  cout<<"a kink"<<endl;
	  break;
	}
	double BzkG;
	int err = TrackModelPhysical::estmateBzkG( mc, mcFwd.x, mcFwd.y, BzkG );
	if( err!=0 ){
	  cout<<"geometry: can not estimate forward Bz. vol "
	      <<h.volume<<" layer "<<h.layer<<" err = "<<err<<endl;
	  //h.Print();
	  //mc.Print();
	  //mcFwd.Print();
	} else {
	  vBz[0][hid] = BzkG;
	}
      } // fwd

      // fit backward field
      {
	if( iBck>=0 && (mc.px*mcBck.px+mc.py*mcBck.py < 0 ) ){
	  cout<<"a kink"<<endl;
	  break;
	}
	double BzkG;
	int err = TrackModelPhysical::estmateBzkG( mc, mcBck.x, mcBck.y, BzkG );
	if( err!=0 ){
	  cout<<"geometry: can not estimate backward Bz, err="<<err<<endl;
	  //h.Print();
	  //mc.Print();
	  //mcBck.Print();
	} else {
	  vBz[2][hid] = BzkG;
	}
      } // backward 
 
      // fit middle field
      if( iFwd>=0 ){	
	TrackModelPhysical tr;
	int err = tr.createXYLast( mcBck.x, mcBck.y, mcBck.z,			      
				   mc.x, mc.y, mc.z,
				   mcFwd.x, mcFwd.y, mcFwd.z, Geo::OriginBzkG );
	if( err!=0 ){
	  cout<<"geometry: can not estimate  Bz for fit, err="<<err<<endl;
	} else {
	  if( tr.pt>0.1 ){
	    vBz[1][hid] = (Geo::OriginBzkG/Geo::CLight ) * mc.pt / tr.pt*mc.q*tr.q;
	  }
	}
      } // middle
    }
  }


  // add measurements to the field fit
  if(1){
    double cut = 3;
    updateLayerField( this, vBz, cut );
  }

  // fill Bz ntuple
  for( unsigned int ipart=0; ipart<mParticles.size(); ipart++ ){          
    Particle &p = mParticles[ipart];             
    for( int ih=0; ih<p.hits.size(); ih++ ){
      int hid = p.hits[ih];
      Hit &h = mHits[hid];
      HitMC &mc = mHitsMC[hid];     
      StatLayer &sl = statLayers[h.layerID];

      //if( vBz[ir][hid]<-100 ) continue;
      float f[100];
      for( int j=0; j<100; j++ ) f[j] = -1.;
      // 0:part :nhits :pt :p :q 5:w :nl :prim :vol :lay 
      //10:layID :mx :my :mz :mr 15:mphi :mt :mpt :mp  :hx 
      //20:hy :hz :hr :hphi :ht 25:Bz0 :Bz1 :Bz2 :fitBz0 :fitBz1
      //30:fitBz2

      f[ 0] = ipart;
      f[ 1] = p.hits.size();
      f[ 2] = p.pt;
      f[ 3] = p.p;
      f[ 4] = p.q;
      f[ 5] = p.w;      
      f[ 6] = p.nLayers;
      f[ 7] = p.prim;
      f[ 8] = h.volume;
      f[ 9] = h.layer;
      f[10] = h.layerID;
      f[11] = mc.x;
      f[12] = mc.y;
      f[13] = mc.z;
      f[14] = mc.r;
      f[15] = mc.phi;
      f[16] = mc.t;
      f[17] = mc.pt;
      f[18] = mc.p;
      f[19] = h.x;
      f[20] = h.y;
      f[21] = h.z;
      f[22] = h.r;
      f[23] = h.phi;
      f[24] = h.t;
      f[25] = vBz[0][hid];
      f[26] = vBz[1][hid];
      f[27] = vBz[2][hid];
      f[28] = sl.getField(0,h.phi, h.t);
      f[29] = sl.getField(1,h.phi, h.t);
      f[30] = sl.getField(2,h.phi, h.t);
      ntGeo->Fill(f);
    }  // particle hits
    //if( ipart%1000==0 ) statFile->Write();
  }// particles


  return;

  // some old stuff 

  for( int iv=0; iv<Geo::NVolumes; iv++){
    Volume &vol = Geo::volumes[iv];    
    for( int il=0; il<vol.nLayers; il++ ){
      int layerId = vol.layerIDs[il];
      //Layer &layer = Geo::layers[layerId];
      FitLayer &flayer = mFitLayers[layerId];

      int nHits = 0;
      double sigmaXY=0, sigmaZ=0;
      int nHitsDense=0;
      double rDense = flayer.mTmin + 0.1*(flayer.mTmax - flayer.mTmin);

      for( int ih=0; ih<mHits.size(); ih++ ){
        Hit &h = mHits[ih];
	HitMC &mc = mHitsMC[ih];
	if( h.layerID!= layerId ) continue;
	if( mc.partID < 0) continue;
	double dx = h.x - mc.x;
	double dy = h.y - mc.y;
	double dz = h.z - mc.z;
	sigmaXY+=dx*dx + dy*dy;
 	sigmaZ+=dz*dz;
	nHits++;
	if( vol.type==0 ){
	  if( fabs(mc.z) <= 0.1*flayer.mTmax ){
	    nHitsDense++;
	  }
	} else {
	  double r = sqrt(mc.x*mc.x+mc.y*mc.y);
	  if( fabs(r) <= rDense ){
	    nHitsDense++;
	  }
	}
      }
      sigmaXY = 3.5*sqrt(sigmaXY/nHits/2);
      sigmaZ = 3.5*sqrt(sigmaZ/nHits);
      double S=0;
      if( vol.type==0 ){	
	S = 2.*TMath::Pi()*flayer.mR*0.1*flayer.mTmax;
      } else {
	S = TMath::Pi()*(rDense*rDense - flayer.mTmin*flayer.mTmin);
      }
      double areaS = S/nHitsDense;
      double areaR = sqrt(areaS/TMath::Pi());
      cout<<"vol "<<iv<<" layer "<<il<<": dens "<< 2*areaR<<" xy=+-"<<sigmaXY<<" z=+-"<<sigmaZ<<endl;
    }
  }
}
