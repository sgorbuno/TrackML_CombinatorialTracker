
/*
  TNtuple *nt = (TNtuple*) gDirectory->FindObjectAny("edgeHits");
  "part:nhits:pt:p:w:side:x:y:z:r"
  TNtuple *nt = (TNtuple*) gDirectory->FindObjectAny("fitHits");
"part:nhits:pt:p:w:fittype:type:x:y:z:r:du:dv:dz:dr:fitpt");
TH1F *h = (TH1F*) gDirectory->FindObjectAny("weights");
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

using namespace std;




Tracker::Tracker()
  :
  mHits(),
  mHitsMC(),
  mParticles()
{

  Geo::init();
  Cuts::init();
}




void Tracker::reconstruct()
{
  cout<<"Reconstruction.."<<endl;
   
  proc.mTracker = this;

  //proc.ReconstructMC();
  proc.vTracklets.clear();
  proc.vTracks.clear();

  
  proc.ReconstructPrim( Cuts::CrCutsPrimV0Pt05, Cuts::PrCutsPrimV0Pt05 );
  proc.ReconstructPrim( Cuts::CrCutsPrimV1Pt05, Cuts::PrCutsPrimV1Pt05 );  
  proc.ReconstructPrim( Cuts::CrCutsPrimV2Pt05, Cuts::PrCutsPrimV1Pt05 );
  
  proc.Reconstruct( Cuts::CrCutsA0Pt05, Cuts::PrCutsA0Pt05 );
  proc.Reconstruct( Cuts::CrCutsB0Pt05, Cuts::PrCutsB0Pt05 );
  //proc.Reconstruct( Cuts::CrCutsD3Pt05, Cuts::PrCutsD3Pt05 );
 
  proc.CleanTracklets( 10 );
  
  proc.ReconstructPrim( Cuts::CrCutsPrimV0Pt03, Cuts::PrCutsPrimV0Pt03 );
  proc.ReconstructPrim( Cuts::CrCutsPrimV1Pt03, Cuts::PrCutsPrimV1Pt03 );  
  proc.ReconstructPrim( Cuts::CrCutsPrimV2Pt03, Cuts::PrCutsPrimV1Pt03 );  
  proc.Reconstruct( Cuts::CrCutsB0Pt03, Cuts::PrCutsB0Pt03 );
  proc.Reconstruct( Cuts::CrCutsD3Pt03, Cuts::PrCutsD3Pt03 );

  proc.CleanTracklets( 10 );
  proc.ReconstructPrim( Cuts::CrCutsPrimV0Pt02, Cuts::PrCutsPrimV0Pt02 );
  proc.ReconstructPrim( Cuts::CrCutsPrimV1Pt02, Cuts::PrCutsPrimV1Pt02 );  
  proc.ReconstructPrim( Cuts::CrCutsPrimV2Pt02, Cuts::PrCutsPrimV1Pt02 );  
 
  proc.CleanTracklets( 8 );
  
  proc.ReconstructPrim( Cuts::CrCutsPrimV0Pt01, Cuts::PrCutsPrimV0Pt01 );
  proc.ReconstructPrim( Cuts::CrCutsPrimV1Pt01, Cuts::PrCutsPrimV1Pt01 );  
  proc.ReconstructPrim( Cuts::CrCutsPrimV2Pt01, Cuts::PrCutsPrimV1Pt01 );  
 

  proc.CleanTracklets( 10 );

  proc.Reconstruct( Cuts::CrCutsPt05Rad, Cuts::PrCutsPt05Rad );

  proc.Reconstruct( Cuts::CrCutsA0Pt03, Cuts::PrCutsA0Pt03 );
  
  proc.Reconstruct( Cuts::CrCutsA5Pt03, Cuts::PrCutsA5Pt03 );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsA5Pt03), Cuts::PrCutsA5Pt03 );

  proc.Reconstruct( Cuts::CrCutsA2Pt01, Cuts::PrCutsA2Pt01 );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsA2Pt01), Cuts::PrCutsA2Pt01 );

  proc.Reconstruct( Cuts::CrCutsB2Pt01, Cuts::PrCutsB2Pt01 );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsB2Pt01), Cuts::PrCutsB2Pt01 );

  proc.Reconstruct( Cuts::CrCutsC2Pt01, Cuts::PrCutsC2Pt01 );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsC2Pt01), Cuts::PrCutsC2Pt01 );

  proc.Reconstruct( Cuts::CrCutsD2Pt01, Cuts::PrCutsD2Pt01 );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsD2Pt01), Cuts::PrCutsD2Pt01 );

  proc.Reconstruct( Cuts::CrCutsE2Pt01, Cuts::PrCutsE2Pt01 );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsE2Pt01), Cuts::PrCutsE2Pt01 );


  proc.Reconstruct( Cuts::CrCutsA3Pt03, Cuts::PrCutsA3Pt03 );
  proc.Reconstruct( Cuts::CrCutsC3Pt03, Cuts::PrCutsC3Pt03 );
  
  proc.Reconstruct( Cuts::CrCutsB3Pt03, Cuts::PrCutsB3Pt03 );

  proc.CleanTracklets( 4 );

  proc.Reconstruct( Cuts::CrCutsB0Pt01, Cuts::PrCutsB0Pt01 );   
  //proc.Reconstruct( Cuts::CrCutsD3Pt01, Cuts::PrCutsD3Pt01 );   
 
  // old stuff which still brings some efficiency 

  proc.Reconstruct( Cuts::CrCutsPt05RadMiddle, Cuts::PrCutsPt05RadMiddle );
  proc.Reconstruct( Cuts::CrCutsPt05FwdMiddle, Cuts::PrCutsPt05FwdMiddle );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsPt05FwdMiddle), Cuts::PrCutsPt05FwdMiddle ); 
 
  
  //proc.CleanTracklets( 7 );  

  proc.Reconstruct( Cuts::CrCutsPt05FwdL, Cuts::PrCutsPt05Fwd );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsPt05FwdL), Cuts::PrCutsPt05Fwd );  


  proc.Reconstruct( Cuts::CrCutsPt05RadMiddleL, Cuts::PrCutsPt05RadMiddle );
  proc.Reconstruct( Cuts::CrCutsPt05FwdMiddleL, Cuts::PrCutsPt05FwdMiddleL );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsPt05FwdMiddleL), Cuts::PrCutsPt05FwdMiddleL ); 
  
  proc.Reconstruct( Cuts::CrCutsPt05FwdC, Cuts::PrCutsPt05Fwd );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsPt05FwdC), Cuts::PrCutsPt05Fwd );

  
  proc.Reconstruct( Cuts::CrCutsPt05FwdMiddleI1, Cuts::PrCutsPt05FwdMiddle );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsPt05FwdMiddleI1), Cuts::PrCutsPt05FwdMiddle );  
  proc.Reconstruct( Cuts::CrCutsPt05MiddleI2, Cuts::PrCutsPt05RadMiddle );
  
  proc.CleanTracklets( 7 );
  
  proc.Reconstruct( Cuts::CrCutsPt01Rad, Cuts::PrCutsPt01Rad );

  proc.CleanTracklets( 7 );  
  
  proc.Reconstruct( Cuts::CrCutsPt01Fwd, Cuts::PrCutsPt01Fwd );
  proc.Reconstruct( Cuts::mirror(Cuts::CrCutsPt01Fwd), Cuts::PrCutsPt01Fwd );
  
  proc.CleanTracklets( 7 );
  //proc.Reconstruct( Cuts::CrCutsA0Pt02, Cuts::PrCutsA0Pt02 );
  //proc.Reconstruct( Cuts::CrCutsA0Pt01, Cuts::PrCutsA0Pt01 );
 
  // proc.Reconstruct( Cuts::CrCutsB0Pt01, Cuts::PrCutsB0Pt01 );
 
  

  /*proc.Reconstruct( Cuts::CrCutsPt01RadMv1, Cuts::PrCutsPt01Rad );
    proc.Reconstruct( Cuts::CrCutsPt01FwdMv1, Cuts::PrCutsPt01Fwd );
    proc.Reconstruct( Cuts::mirror(Cuts::CrCutsPt01FwdMv1), Cuts::PrCutsPt01Fwd );
  */

  //proc.CleanTracklets( 3 );
  //proc.vTracklets.clear();  


  //proc.TrackletFitTestKF( Cuts::CrCutsA5Pt03 );
 

  //proc.Reconstruct( Cuts::CrCutsB0Pt05, Cuts::PrCutsB0Pt05 );
  //proc.CleanTracklets( 1 );
  
  //proc.Reconstruct( Cuts::CrCutsA5Pt03, Cuts::PrCutsA5Pt03 );
 

  proc.CleanTracklets( 1 );
  proc.TrackletEfficiency( 0, -0.3, -0.3);
  
  //analyzeGeometry();
}
