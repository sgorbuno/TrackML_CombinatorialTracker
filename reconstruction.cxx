
/*
  TNtuple *nt = (TNtuple*) gDirectory->FindObjectAny("edgeHits");
  "part:nhits:pt:p:w:side:x:y:z:r"
  TNtuple *nt = (TNtuple*) gDirectory->FindObjectAny("fitHits");
"part:nhits:pt:p:w:fittype:type:x:y:z:r:du:dv:dz:dr:fitpt");
TH1F *h = (TH1F*) gDirectory->FindObjectAny("weights");
 */

#include "Tracker.h"
#include "TString.h"
#include "TFile.h"
#include <iostream>

using namespace std;

int main()
{  
  bool analyseTruth = false;

  const int nEvents = 1;
  const int firstEvent=1000;
  //TString dir = "../task/train_100_events/";  
  TString dir = "data/train_100_events/";  
  
  /*
  const int nEvents = 125;
  const int firstEvent=0;
  TString dir = "data/test/";  
  */

  /*
  const int nEvents = 1000;
  const int firstEvent=1000;
  TString dir = "../task/train_1/";  
  */

  ofstream out("mysubmission.csv");    
  if( !out.is_open() ){
    cout<<"Can not open output file"<<endl;
    exit(0);
  }
  
  out<<"event_id,hit_id,track_id"<<endl;
  
  Tracker tracker;
  
  long int currentID=1;

  for( int event = firstEvent; event<firstEvent+nEvents; event++){    
    cout<<"read event "<<event<<endl;
    tracker.readEvent( dir.Data(),  event, analyseTruth );
    
    //tracker.analyzeGeometry(0);    
    //continue;

    tracker.reconstruct();

    cout<<"Event "<<event<<": reconstructed "<<tracker.proc.vTracks.size()<<" tracks"<<endl;

    TrackletConstructor &proc = tracker.proc;

    for( int ih=0; ih<tracker.mHits.size(); ih++ ){
      Hit &h = tracker.mHits[ih];
      h.trackID=0;
    }

    currentID=1;

    for( int itr=0; itr<proc.vTracks.size(); itr++ ){
      Tracklet &t = proc.vTracks[itr];
      if( t.vHits.size()==0 ) continue;
      for( int ih=0; ih<t.vHits.size(); ih++ ){
	Hit &h = tracker.mHits[t.vHits[ih]];
	if( h.trackID>0 ){
	  cout<<"reconstruction.cxx: Wrong hit indexing!"<<endl;
	  cout<<"track "<<itr<<" nhits "<<t.vHits.size()<<endl;
	  cout<<"hit "<<ih<<" trackID "<<h.trackID<<endl;
	  h.Print();
	  exit(1);
	}
	h.trackID=currentID;
      }
      currentID++;
    }
    
    for( int ih=0; ih<tracker.mHits.size(); ih++ ){
      Hit &h = tracker.mHits[ih];
      out<<event<<","<<ih+1<<","<<h.trackID<<endl;
    }
  } // events
  
  //tracker.analyzeGeometry(1);

  out.close();
  return 0;
}
