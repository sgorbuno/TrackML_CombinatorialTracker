#include "Cuts.h"

// === primary 

CreationPrimCuts Cuts::CrCutsPrimV0Pt05;
ProlongationCuts Cuts::PrCutsPrimV0Pt05;
CreationPrimCuts Cuts::CrCutsPrimV0Pt03;
ProlongationCuts Cuts::PrCutsPrimV0Pt03;
CreationPrimCuts Cuts::CrCutsPrimV0Pt02;
ProlongationCuts Cuts::PrCutsPrimV0Pt02;
CreationPrimCuts Cuts::CrCutsPrimV0Pt01;
ProlongationCuts Cuts::PrCutsPrimV0Pt01;

CreationPrimCuts Cuts::CrCutsPrimV1Pt05;
CreationPrimCuts Cuts::CrCutsPrimV2Pt05;
ProlongationCuts Cuts::PrCutsPrimV1Pt05;

CreationPrimCuts Cuts::CrCutsPrimV1Pt03;
CreationPrimCuts Cuts::CrCutsPrimV2Pt03;
ProlongationCuts Cuts::PrCutsPrimV1Pt03;

CreationPrimCuts Cuts::CrCutsPrimV1Pt02;
CreationPrimCuts Cuts::CrCutsPrimV2Pt02;
ProlongationCuts Cuts::PrCutsPrimV1Pt02;

CreationPrimCuts Cuts::CrCutsPrimV1Pt01;
CreationPrimCuts Cuts::CrCutsPrimV2Pt01;
ProlongationCuts Cuts::PrCutsPrimV1Pt01;

CreationCuts     Cuts::CrCutsA0Pt05;
ProlongationCuts Cuts::PrCutsA0Pt05;
CreationCuts     Cuts::CrCutsA0Pt03;
ProlongationCuts Cuts::PrCutsA0Pt03;
CreationCuts     Cuts::CrCutsA0Pt02;
ProlongationCuts Cuts::PrCutsA0Pt02;
CreationCuts     Cuts::CrCutsA0Pt01;
ProlongationCuts Cuts::PrCutsA0Pt01;

CreationCuts     Cuts::CrCutsB0Pt05;
ProlongationCuts Cuts::PrCutsB0Pt05;
CreationCuts     Cuts::CrCutsB0Pt03;
ProlongationCuts Cuts::PrCutsB0Pt03;
CreationCuts     Cuts::CrCutsB0Pt01;
ProlongationCuts Cuts::PrCutsB0Pt01;

CreationCuts     Cuts::CrCutsA2Pt01;
ProlongationCuts Cuts::PrCutsA2Pt01;
CreationCuts     Cuts::CrCutsB2Pt01;
ProlongationCuts Cuts::PrCutsB2Pt01;
CreationCuts     Cuts::CrCutsC2Pt01;
ProlongationCuts Cuts::PrCutsC2Pt01;
CreationCuts     Cuts::CrCutsD2Pt01;
ProlongationCuts Cuts::PrCutsD2Pt01;

CreationCuts     Cuts::CrCutsE2Pt01;
ProlongationCuts Cuts::PrCutsE2Pt01;

CreationCuts     Cuts::CrCutsA3Pt03;
ProlongationCuts Cuts::PrCutsA3Pt03;
CreationCuts     Cuts::CrCutsC3Pt03;
ProlongationCuts Cuts::PrCutsC3Pt03;

CreationCuts     Cuts::CrCutsB3Pt03;
ProlongationCuts Cuts::PrCutsB3Pt03;

CreationCuts     Cuts::CrCutsD3Pt05;
ProlongationCuts Cuts::PrCutsD3Pt05;
CreationCuts     Cuts::CrCutsD3Pt03;
ProlongationCuts Cuts::PrCutsD3Pt03;
CreationCuts     Cuts::CrCutsD3Pt01;
ProlongationCuts Cuts::PrCutsD3Pt01;


CreationCuts     Cuts::CrCutsA5Pt03;
ProlongationCuts Cuts::PrCutsA5Pt03;




  //==========   Radial cuts
  
CreationCuts Cuts::CrCutsPt1Rad;
CreationCuts Cuts::CrCutsPt05Rad;
CreationCuts Cuts::CrCutsPt05RadMiddle;
CreationCuts Cuts::CrCutsPt05RadMiddleL;
CreationCuts Cuts::CrCutsPt05MiddleI2;
CreationCuts Cuts::CrCutsPt03Rad;
CreationCuts Cuts::CrCutsPt01Rad;
CreationCuts Cuts::CrCutsPt01MiddleRad;
CreationCuts Cuts::CrCutsPtOld11MiddleRad;
CreationCuts Cuts::CrCutsPt01RadMv1;
 

  //==========  Forward cuts

CreationCuts Cuts::CrCutsPt1Fwd; 
CreationCuts Cuts::CrCutsPt05Fwd;
CreationCuts Cuts::CrCutsPt05FwdL;
CreationCuts Cuts::CrCutsPt05FwdC;
CreationCuts Cuts::CrCutsPt05FwdMiddle;
CreationCuts Cuts::CrCutsPt05FwdMiddleL;
CreationCuts Cuts::CrCutsPt05FwdMiddleI1;
CreationCuts Cuts::CrCutsPt03Fwd;
CreationCuts Cuts::CrCutsPt01Fwd;
CreationCuts Cuts::CrCutsPt1MiddleFwd;
CreationCuts Cuts::CrCutsPt01MiddleFwd;

  // ===========   prolongation cuts ====
  
ProlongationCuts Cuts::PrCutsPt1Fwd ;
ProlongationCuts Cuts::PrCutsPt05Fwd ;
ProlongationCuts Cuts::PrCutsPt05FwdMiddle ;
ProlongationCuts Cuts::PrCutsPt05FwdMiddleL ;
ProlongationCuts Cuts::PrCutsPt05FwdMiddleI1 ;
ProlongationCuts Cuts::PrCutsPt01Fwd ;

ProlongationCuts Cuts::PrCutsPt05RadLoose ;
ProlongationCuts Cuts::PrCutsPt05Rad ;
ProlongationCuts Cuts::PrCutsPt05RadMiddle ;
//ProlongationCuts Cuts::PrCutsPt05I1;
ProlongationCuts Cuts::PrCutsPt03Rad;
ProlongationCuts Cuts::PrCutsPt01Rad ; 
ProlongationCuts Cuts::PrCutsPt01MiddleRad ;


CreationCuts Cuts::mirror( const CreationCuts  &fwd ) 
{
  CreationCuts bck = fwd;
  for( int i=0; i<3; i++){
    const Layer &l = Geo::layers[fwd.layerCuts[i].iLayer];
    int vmirror = l.volume;
    switch( l.volume ){
    case 1:
    case 4:
    case 7:
      vmirror = l.volume+1;
      break;
    case 2:
    case 5:
    case 8:
      vmirror = l.volume-1;
      break;
    case 0:
    case 3:
    case 6:
      vmirror = l.volume;
      break;
    default:
      cout<<"cuts mirroring: something wrong!!"<<endl;
      exit(1);
    }    
    bck.layerCuts[i].iLayer = Geo::getLayerID(vmirror,l.layer);
  }
  return bck;
}



void Cuts::init()
{
  
  //  ==== Primary cuts
    
  CrCutsPrimV0Pt05 = 
    {
      .43, // pt
      {  //    layer               phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),  0.03,  1.6,  1.6,     -1,   1.3 }
      }
    };
  
  PrCutsPrimV0Pt05 = { 
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.005, 0.1,  0.06, 0.1,  0.0050, 0.0045 }, // vol 0
      { 0.005, 0.03, 0.06, 0.03, 0.0050, 0.0012 }, // vol 1,2
      { 0.005, 0.03, 0.06, 0.03, 0.0050, 0.0012 },  
      { 0.01,  0.4,  0.23,  0.35, 0.01,   0.11   }, // vol 3
      { 0.01,  0.3,  0.4,  0.3,   0.03,   0.11   }, // vol 4,5
      { 0.01,  0.3,  0.4,  0.3,   0.03,   0.11   }, 
      { 0.007, 1.9,  0.3,  1.8,   0.016,  1.1    }, // vol 6
      { 0.03,  1.6,   2.,  1.7,   0.3,    0.9    }, // vol 7,8
      { 0.03,  1.6,   2.,  1.7,   0.3,    0.9    }
    }
  };
 

 CrCutsPrimV0Pt03 = 
    {
      .25, // pt
      {  //    layer               phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),  0.05,  1.6,  1.6,     -1,   1.3 }
      }
    };
  
  PrCutsPrimV0Pt03 = { 
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.009, 0.16,  0.11, 0.16,  0.005,  0.02 }, // vol 0
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 },  
      { 0.025,  0.8,  0.81, 0.6,   0.02,  0.11   }, // vol 3
      { 0.03,  0.8,  0.7,  0.75,   0.04,    0.11   }, // vol 4,5
      { 0.03,  0.8,  0.7,  0.75,   0.04,    0.11   }, 
      { 0.035, 1.4,  2.6,    1.6, 0.1,   1.2    }, // vol 6
      { 0.025,  2.5,   1.9,   2.5,   0.9,    0.9    }, // vol 7,8
      { 0.025,  2.5,   1.9,   2.5,   0.9,    0.9    }
    }
  };
  
  CrCutsPrimV0Pt02 = 
    {
      .17, // pt
      {  //    layer               phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),  0.07,  1.7,  1.7,     -1,   1.4 }
      }
    };
  
  PrCutsPrimV0Pt02 = { 
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.014, 0.22,  0.20, 0.22,  0.0045,  0.05 }, // vol 0
      { 0.009, 0.12,  0.15, 0.12,  0.006, 0.015 }, // vol 1,2
      { 0.009, 0.12,  0.15, 0.12,  0.006, 0.015 },
      { 0.035,  1.6,  1.3,  1.3,   0.06,  0.11   }, // vol 3
      { 0.019,  0.9,  0.9,   1.,   0.07,   0.11   }, // vol 4,5
      { 0.019,  0.9,  0.9,   1.,   0.07,   0.11   }, 
      { 0.04,  1.6,  2.6,    1.6,  0.06,   1.2    }, // vol 6
      { 0.2,  4.0,   2.5,   4.,   1.,    0.9    }, // vol 7,8
      { 0.2,  4.0,   2.5,   4.,   1.,    0.9    }
    }
  };

  CrCutsPrimV0Pt01 = 
    {
      .09, // pt
      {  //    layer               phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),  0.095,  1.7,  1.7,     -1,   1.4 }
      }
    };
  
  PrCutsPrimV0Pt01 = { 
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.015, 0.35,  0.21, 0.35,  0.0065,  0.05 }, // vol 0
      { 0.012, 0.12,  0.16, 0.12,  0.006, 0.022 }, // vol 1,2
      { 0.012, 0.12,  0.16, 0.12,  0.006, 0.022 },
      { 0.05,  2.,    1.3,  1.3,   0.06,  0.12   }, // vol 3
      { 0.06,  1.5,  2.1,   1.2,   0.08,   0.11   }, // vol 4,5
      { 0.06,  1.5,  2.1,   1.2,   0.08,   0.11   }, 
      { 0.04,  1.6,  2.6,    1.6,  0.06,   1.2    }, // vol 6
      { 0.2,  4.0,   2.5,   4.,   1.,    0.9    }, // vol 7,8
      { 0.2,  4.0,   2.5,   4.,   1.,    0.9    }
    }
  };


  CrCutsPrimV1Pt05 = 
    {
      .3, // pt
      {  //    layer               phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 1, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 1, 1),  0.015, 0.06,  0.06,   -1,   1.6 }
      }
    };

  CrCutsPrimV2Pt05 = CrCutsPrimV1Pt05;
    {
      CrCutsPrimV2Pt05.layerCuts[0].iLayer = Geo::getLayerID( 2, 0);
      CrCutsPrimV2Pt05.layerCuts[1].iLayer = Geo::getLayerID( 2, 1);
    }
  
  PrCutsPrimV1Pt05 = { 
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.005, 0.31, 0.016, 0.31,  0.0035,0.08 }, // vol 0
      { 0.005, 0.03, 0.06, 0.03, 0.005, 0.0045 }, // vol 1,2
      { 0.005, 0.03, 0.06, 0.03, 0.005, 0.0045 },  
      { 0.01,  0.4,  0.23,  0.35, 0.01,   0.11   }, // vol 3
      { 0.01,  0.3,  0.4,  0.3,   0.03,   0.11   }, // vol 4,5
      { 0.01,  0.3,  0.4,  0.3,   0.03,   0.11   }, 
      { 0.007, 1.9,  0.3,  1.8,   0.016,  1.1    }, // vol 6
      { 0.03,  1.6,   2.,  1.7,   0.3,    0.9    }, // vol 7,8
      { 0.03,  1.6,   2.,  1.7,   0.3,    0.9    }
    }
  };


 
  CrCutsPrimV1Pt03 = 
    {
      .2, // pt
      {  //    layer               phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 1, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 1, 1),  0.022, 0.1,  0.1,   -1,   3 } // main duz 1.6
      }
    };

  CrCutsPrimV2Pt03 = CrCutsPrimV1Pt03;
  {
    CrCutsPrimV2Pt03.layerCuts[0].iLayer = Geo::getLayerID( 2, 0);
    CrCutsPrimV2Pt03.layerCuts[1].iLayer = Geo::getLayerID( 2, 1);
  }
  
  PrCutsPrimV1Pt03 = { 
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.01, 0.7,   0.03,  0.7,  0.0035,0.08 }, // vol 0
      { 0.005, 0.05, 0.06, 0.05, 0.005, 0.0045 }, // vol 1,2
      { 0.005, 0.05, 0.06, 0.05, 0.005, 0.0045 },  
      { 0.01,  0.4,  0.23,  0.35, 0.01,   0.11   }, // vol 3
      { 0.025,  0.8,  0.6,  0.8,  0.03,   0.11   }, // vol 4,5
      { 0.025,  0.8,  0.6,  0.8,  0.03,   0.11   }, 
      { 0.007, 1.9,  0.3,  1.8,   0.016,  1.1    }, // vol 6
      { 0.03,  1.6,   2.,  1.7,   0.3,    0.9    }, // vol 7,8
      { 0.03,  1.6,   2.,  1.7,   0.3,    0.9    }
    }
  };


   CrCutsPrimV1Pt02 = 
    {
      0.14, // pt
      {  //    layer               phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 1, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 1, 1),  0.04,  0.1,  0.1,   -1,   3 } // main duz 1.6
      }
    };

  CrCutsPrimV2Pt02 = CrCutsPrimV1Pt02;
  {
    CrCutsPrimV2Pt02.layerCuts[0].iLayer = Geo::getLayerID( 2, 0);
    CrCutsPrimV2Pt02.layerCuts[1].iLayer = Geo::getLayerID( 2, 1);
  }
  
  PrCutsPrimV1Pt02 = { 
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.012, 0.7,   0.04,  0.7,  0.0035,0.09 }, // vol 0
      { 0.005, 0.05, 0.06, 0.05, 0.005, 0.0045 }, // vol 1,2
      { 0.005, 0.05, 0.06, 0.05, 0.005, 0.0045 },  
      { 0.01,  0.4,  0.23,  0.35, 0.01,   0.11   }, // vol 3
      { 0.025,  0.8,  0.4,  0.65,   0.03,   0.11   }, // vol 4,5
      { 0.025,  0.8,  0.4,  0.65,   0.03,   0.11   }, 
      { 0.007, 1.9,  0.3,  1.8,   0.016,  1.1    }, // vol 6
      { 0.03,  1.6,   2.,  1.7,   0.3,    0.9    }, // vol 7,8
      { 0.03,  1.6,   2.,  1.7,   0.3,    0.9    }
    }
  };

   CrCutsPrimV1Pt01 = 
    {
      -0.14, // pt
      {  //    layer               phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 1, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 1, 1),  0.05,  0.13,  0.13,   -1,   3 } // main duz 1.6
      }
    };

  CrCutsPrimV2Pt01 = CrCutsPrimV1Pt01;
  {
    CrCutsPrimV2Pt01.layerCuts[0].iLayer = Geo::getLayerID( 2, 0);
    CrCutsPrimV2Pt01.layerCuts[1].iLayer = Geo::getLayerID( 2, 1);
  }
  
  PrCutsPrimV1Pt01 = { 
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.012, 0.7,   0.04,  0.7,  0.0035,0.09 }, // vol 0
      { 0.0065, 0.09, 0.08, 0.08, 0.005, 0.0045 }, // vol 1,2
      { 0.0065, 0.09, 0.08, 0.08, 0.005, 0.0045 },  
      { 0.01,  0.4,  0.23,  0.35, 0.01,   0.11   }, // vol 3
      { 0.03,  0.9,  0.8,  0.9,   0.04,   0.12   }, // vol 4,5
      { 0.03,  0.9,  0.8,  0.9,   0.04,   0.12   }, 
      { 0.007, 1.9,  0.3,  1.8,   0.016,  1.1    }, // vol 6
      { 0.03,  1.6,   2.,  1.7,   0.3,    0.9    }, // vol 7,8
      { 0.03,  1.6,   2.,  1.7,   0.3,    0.9    }
    }
  };

 
  // ====== A0 ===
  /*
  CrCutsA0Pt05 = 
    {
      0, // use vertex
      -.4, // pt
      {  //    layer               phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),   0.2,  18.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 0, 2),  0.07,  0.1,   -1,     -1,   0.07 }
      }
    };
  */
  
  CrCutsA0Pt05 = 
    {
      1, // use vertex
      .4, // pt
      {  //    layer               phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),   0.2,  3.6,   3.6,     -1,    -1 },
 	{ Geo::getLayerID( 0, 2),  0.021,  0.1,  0.15,     -1,   0.07 }
      }
    };
  
  PrCutsA0Pt05 = { //PrCutsPrimV0Pt05;//PrCutsA0Pt03;
     {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.015, 0.12,  0.05, 0.06,  0.005,  0.03 }, // vol 0
      { 0.009, 0.04,  0.05, 0.035, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.04,  0.05, 0.035, 0.0045, 0.015 },  
      { 0.008,  0.3,  0.2,  0.3,   0.01,   0.11   }, // vol 3
      { 0.025,  0.3,  0.2,  0.3,   0.04,    0.11   }, // vol 4,5
      { 0.025,  0.3,  0.2,  0.3,   0.04,    0.11   }, 
      { 0.006,  1.5,  0.4,   1.5,   0.01,   1.2    }, // vol 6
      { 0.01,   1.5,  0.8,   1.5,   0.3,    0.9    }, // vol 7,8
      { 0.01,   1.5,  0.8,   1.5,   0.3,    0.9    }
    }
  };

  
  CrCutsA0Pt03 = 
    {
      0, // use vertex
      .25, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),   0.5,  15.5,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 0, 2),   0.2,   0.5,   -1,     -1,   0.15 }
      }
    };
  
  PrCutsA0Pt03 = {
   {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.022, 0.7,  0.11, 0.16,  0.005,  0.02 }, // vol 0
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 },  
      { 0.025,  0.8,  0.5, 0.6,   0.02,  0.11   }, // vol 3
      { 0.03,  0.8,  0.7,  0.75,   0.04,    0.11   }, // vol 4,5
      { 0.03,  0.8,  0.7,  0.75,   0.04,    0.11   }, 
      { 0.035, 1.4,  2.6,    1.6, 0.1,   1.2    }, // vol 6
      { 0.025,  2.5,   1.9,   2.5,   0.9,    0.9    }, // vol 7,8
      { 0.025,  2.5,   1.9,   2.5,   0.9,    0.9    }
    }
  };


 
  CrCutsA0Pt02 = 
    {
      0, // use vertex
      .15, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),   1.1,  15.5,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 0, 2),   0.2,   0.5,   -1,     -1,   0.3 }
      }
    };

  PrCutsA0Pt02 = PrCutsPrimV0Pt02;   
 
  CrCutsA0Pt01 = 
    {
      0, // use vertex
      -.15, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),   0.6,  15.5,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 0, 2),   0.2,   0.55,   -1,     -1,   0.15 }
      }
    };

  PrCutsA0Pt01 = PrCutsPrimV0Pt01;   


  //==== B0

  
  CrCutsB0Pt05 = 
    {
      0, // use vertex
      .48, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 0, 2),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 3),   0.25,    17.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 3, 0),   0.15,   0.25,   -1,     -1,   0.1 }
      }
      };

  PrCutsB0Pt05 = {//PrCutsPrimV0Pt03;   
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.012, 0.4,  0.04, 0.08,  0.005,  0.02 }, // vol 0
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 },  
      { 0.006,  0.3,  0.2, 0.4,   0.01,  0.11   }, // vol 3
      { 0.003,  0.2,  0.2,  0.2,   0.04,    0.11   }, // vol 4,5
      { 0.003,  0.2,  0.2,  0.2,   0.04,    0.11   }, 
      { 0.007,  1.5,  0.4,    1.5, 0.01,   1.2    }, // vol 6
      { 0.005,  1.7,   0.6,   1.8,   0.3,    0.9    }, // vol 7,8
      { 0.005,  1.7,   0.6,   1.8,   0.3,    0.9    }
    }
  };


 CrCutsB0Pt03 = 
    {
      0, // use vertex
      .14, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 0, 2),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 3),   0.3,    17.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 3, 0),   0.25,   0.65,   -1,     -1,   0.16 }
      }
    };

  PrCutsB0Pt03 = {//PrCutsPrimV0Pt03;   
    {  // phi   T      V     UZ    dupV    dupUZ    
      { 0.08,  0.7,  0.07, 0.16,  0.005,  0.02 }, // vol 0
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 },  
      { 0.012,  0.8,  0.5, 0.45,   0.01,  0.11   }, // vol 3
      { 0.08,   0.8,  0.5,  0.6,   0.01,    0.11   }, // vol 4,5
      { 0.08,   0.8,  0.5,  0.6,   0.01,    0.11   }, 
      { 0.03,   2.3,  1.5,    2.2,  0.1,   1.2    }, // vol 6
      { 0.016,  2.5,   1.9,   2.8,   0.3,    0.9    }, // vol 7,8
      { 0.016,  2.5,   1.9,   2.8,   0.3,    0.9    }
    }
  };

  CrCutsB0Pt01 = 
    {
      0, // use vertex
      .14, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 0, 2),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 3),   0.45,    17.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 3, 0),   0.25,   1.8,   -1,     -1,   0.2 }
      }
    };

  PrCutsB0Pt01 = {
    {  // phi   T      V     UZ    dupV    dupUZ    
      { 0.08,  0.7,  0.01, 0.2,  0.005,  0.02 }, // vol 0
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 },  
      { 0.08,  2.0,   2., 0.5,   0.01,  0.11   }, // vol 3
      { 0.08,   2.5,  2.5,  1.,   0.01,    0.11   }, // vol 4,5
      { 0.08,   2.5,  2.5,  1.,   0.01,    0.11   }, 
      { 0.08,   2.3,  2.,    2.2,  0.1,   1.2    }, // vol 6
      { 0.04,  2.5,   3.,   2.8,   0.3,    0.9    }, // vol 7,8
      { 0.04,  2.5,   3.,   2.8,   0.3,    0.9    }
    }
  };

  
  // A2 

  CrCutsA2Pt01 = 
    {
      1, // use vertex
      -.1, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 2, 2),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 3),   0.4,  1.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 2, 4),   0.2,  0.1,   -1,     -1,   0.04 }
      }
    };
  
  PrCutsA2Pt01 = PrCutsPrimV1Pt01;   
 
  // B2 
  
  CrCutsB2Pt01 = 
    {
      1, // use vertex
      -.1, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 2, 3),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 4),   0.4,  1.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 2, 5),   0.2,  0.1,   -1,     -1,   0.04 }
      }
    };
  
  PrCutsB2Pt01 = PrCutsA2Pt01;
 
 // C2 
  
  CrCutsC2Pt01 = 
    {
      1, // use vertex
      -.1, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 2, 4),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 5),   0.4,  1.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 2, 6),   0.2,  0.1,   -1,     -1,   0.04 }
      }
    };
  
  PrCutsC2Pt01 = PrCutsA2Pt01;
 
 CrCutsD2Pt01 = 
    {
      1, // use vertex
      -.1, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 2, 1),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 2),   0.4,  1.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 2, 3),   0.2,  0.1,   -1,     -1,   0.04 }
      }
    };
  
  PrCutsD2Pt01 = PrCutsA2Pt01;

  CrCutsE2Pt01 = 
    {
      1, // use vertex
      -.1, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 2, 0),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 1),   0.4,  1.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 2, 2),   0.2,  0.1,   -1,     -1,   0.04 }
      }
    };
  
  PrCutsE2Pt01 = PrCutsA2Pt01;

  //==== A3


  CrCutsA3Pt03 = 
    {
      0, // use vertex
      .25, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 3, 0),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 1),   0.5,  15.5,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 3, 2),   0.06,   0.5,   -1,     -1,   0.25 }
      }
    };
  
  PrCutsA3Pt03 = {//PrCutsPrimV0Pt03;   
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.022, 0.7,  0.11, 0.16,  0.005,  0.02 }, // vol 0
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 },  
      { 0.025,  0.8,  0.5, 0.6,   0.02,  0.11   }, // vol 3
      { 0.03,  0.8,  0.7,  0.75,   0.04,    0.11   }, // vol 4,5
      { 0.03,  0.8,  0.7,  0.75,   0.04,    0.11   }, 
      { 0.035, 1.4,  2.6,    1.6, 0.1,   1.2    }, // vol 6
      { 0.025,  2.5,   1.9,   2.5,   0.9,    0.9    }, // vol 7,8
      { 0.025,  2.5,   1.9,   2.5,   0.9,    0.9    }
    }
  };

 //==== C3


  CrCutsC3Pt03 = 
    {
      0, // use vertex
      .25, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 3, 1),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 2),   0.6,  30.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 3, 3),   0.22,   1.2,   -1,     -1,   0.25 }
      }
    };
  
  PrCutsC3Pt03 = PrCutsA3Pt03;


  //==== B3

  
 CrCutsB3Pt03 = 
    {
      1, // use vertex
      .25, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 0, 3),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 0),   0.4,   5.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 3, 1),   0.3,   1.2,   -1,     -1,   0.25 }
      }
    };

  PrCutsB3Pt03 = {//PrCutsPrimV0Pt03;   
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.022, 0.7,  0.11, 0.16,  0.005,  0.02 }, // vol 0
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 },  
      { 0.025,  0.8,  0.5, 0.6,   0.02,  0.11   }, // vol 3
      { 0.03,  0.8,  0.7,  0.75,   0.04,    0.11   }, // vol 4,5
      { 0.03,  0.8,  0.7,  0.75,   0.04,    0.11   }, 
      { 0.035, 1.4,  2.6,    1.6, 0.1,   1.2    }, // vol 6
      { 0.025,  2.5,   1.9,   2.5,   0.9,    0.9    }, // vol 7,8
      { 0.025,  2.5,   1.9,   2.5,   0.9,    0.9    }
    }
  };



  //==== D3

  
  CrCutsD3Pt05 = 
    {
      0, // use vertex
      .48, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 0, 3),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 0),   0.25,    17.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 3, 1),   0.15,   0.25,   -1,     -1,   0.1 }
      }
      };

  PrCutsD3Pt05 = {//PrCutsPrimV0Pt03;   
    {  // phi  T       V     UZ    dupV    dupUZ    
      { 0.012, 0.4,  0.04, 0.08,  0.005,  0.02 }, // vol 0
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 },  
      { 0.006,  0.3,  0.2, 0.4,   0.01,  0.11   }, // vol 3
      { 0.003,  0.2,  0.2,  0.2,   0.04,    0.11   }, // vol 4,5
      { 0.003,  0.2,  0.2,  0.2,   0.04,    0.11   }, 
      { 0.007,  1.5,  0.4,    1.5, 0.01,   1.2    }, // vol 6
      { 0.005,  1.7,   0.6,   1.8,   0.3,    0.9    }, // vol 7,8
      { 0.005,  1.7,   0.6,   1.8,   0.3,    0.9    }
    }
  };


 CrCutsD3Pt03 = 
    {
      0, // use vertex
      .14, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 0, 3),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 0),   0.3,    17.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 3, 1),   0.25,   0.65,   -1,     -1,   0.16 }
      }
    };

  PrCutsD3Pt03 = {//PrCutsPrimV0Pt03;   
    {  // phi   T      V     UZ    dupV    dupUZ    
      { 0.08,  0.7,  0.07, 0.16,  0.005,  0.02 }, // vol 0
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 },  
      { 0.012,  0.8,  0.5, 0.45,   0.01,  0.11   }, // vol 3
      { 0.08,   0.8,  0.5,  0.6,   0.01,    0.11   }, // vol 4,5
      { 0.08,   0.8,  0.5,  0.6,   0.01,    0.11   }, 
      { 0.03,   2.3,  1.5,    2.2,  0.1,   1.2    }, // vol 6
      { 0.016,  2.5,   1.9,   2.8,   0.3,    0.9    }, // vol 7,8
      { 0.016,  2.5,   1.9,   2.8,   0.3,    0.9    }
    }
  };

  CrCutsD3Pt01 = 
    {
      0, // use vertex
      .14, // pt
      {  //    layer               phiL     TL    TH      dv    duz
	{ Geo::getLayerID( 0, 3),    -1,    -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 0),   0.45,    17.,   -1,     -1,    -1 },
 	{ Geo::getLayerID( 3, 1),   0.25,   1.8,   -1,     -1,   0.2 }
      }
    };

  PrCutsD3Pt01 = {
    {  // phi   T      V     UZ    dupV    dupUZ    
      { 0.08,  0.7,  0.01, 0.2,  0.005,  0.02 }, // vol 0
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 }, // vol 1,2
      { 0.009, 0.05,  0.1, 0.045, 0.0045, 0.015 },  
      { 0.08,  2.0,   2., 0.5,   0.01,  0.11   }, // vol 3
      { 0.08,   2.5,  2.5,  1.,   0.01,    0.11   }, // vol 4,5
      { 0.08,   2.5,  2.5,  1.,   0.01,    0.11   }, 
      { 0.08,   2.3,  2.,    2.2,  0.1,   1.2    }, // vol 6
      { 0.04,  2.5,   3.,   2.8,   0.3,    0.9    }, // vol 7,8
      { 0.04,  2.5,   3.,   2.8,   0.3,    0.9    }
    }
  };




  // A5


  CrCutsA5Pt03 =  
    {
      1, // use vertex
      .1, // pt
      { //    layer                phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 3, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 5, 0),  0.4,   1.5,   -1,     -1,    -1 },
	{ Geo::getLayerID( 5, 1),  0.06,  1.2,  - 1,     -1,  0.7 }
      }
    };


  PrCutsA5Pt03 = { 
    {  // phi    T       V     UZ     dupV   dupUZ
      { 0.015,  0.5,   0.2,   0.6,  0.005,  0.025  }, // vol 0
      { 0.02,  0.22,   0.3,  0.35,  0.005,  0.002  }, // vol 1,2
      { 0.02,  0.22,   0.3,  0.35,  0.005,  0.002  },  
      { 0.01,   1.0,   0.6,  1.2,  0.008,  0.08  }, // vol 3
      { 0.008,  0.41,  0.3,   0.5,   0.04,    0.1  }, // vol 4,5
      { 0.008,  0.41,  0.3,   0.5,   0.04,    0.1  }, 
      { 0.005,  2.0,   0.6,    3.,    0.5,    5.0  }, // vol 6 - empty
      { 0.03,   2.,   1.5,   3.,    0.5,    1.  }, // vol 7,8 
      { 0.03,   2.,   1.5,   3.,    0.5,    1.  }
    }
  };



 //==========   Radial cuts
  
  CrCutsPt1Rad = 
    {
       1, // use vertex
     0.9, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),   0.1,  2.0,  2.0,     -1,    -1 },
	{ Geo::getLayerID( 0, 2),   0.1,  1.6,  1.6,  0.001,  0.05 }
      }
    };

  CrCutsPt05Rad = 
    {
       1, // use vertex
     .45, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),   0.1,  2.1,  2.0,     -1,    -1 },
	{ Geo::getLayerID( 0, 2),  0.023,  0.6,  0.15,  0.001,  0.06 }
      }
    };

  CrCutsPt05RadMiddle = 
    {
       1, // use vertex
     .45, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 3, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 1),  0.08,   3.,  3.0,     -1,    -1 },
	{ Geo::getLayerID( 3, 2),   0.1,  1.2,  0.4,  0.001,  0.1 }
      }
    };

  CrCutsPt05RadMiddleL = 
    {
       1, // use vertex
     .45, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 3, 1),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 2),  0.08,   3.,  3.0,     -1,    -1 },
	{ Geo::getLayerID( 3, 3),   0.1,  1.2,  0.4,  0.001,  0.1 }
      }
    };

 CrCutsPt05MiddleI2 = 
    {
       1, // use vertex
     .45, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 3, 2),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 3),  0.12,   3.,  2.2,     -1,    -1 },
	{ Geo::getLayerID( 6, 0),  0.12,  1.2,  1.2,  0.001,  0.32 }
      }
    };


  CrCutsPt03Rad = 
    {
       1, // use vertex
     0.25, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),   0.4,  20.,  20.,     -1,    -1 },
	{ Geo::getLayerID( 0, 2),   0.04,  0.62,  0.2,  0.001,    0.08 }
      }
    };

  CrCutsPt01Rad = 
    {
       1, // use vertex
     -1., // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 0, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 0, 1),   .4,  20.,  20.,     -1,    -1 },
	{ Geo::getLayerID( 0, 2),   0.08,  0.8,  0.6,  0.001,   0.1 }
      }
    };
  


  CrCutsPt01MiddleRad = 
    {
       1, // use vertex
     0.5, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 3, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 1),    1.,  80.0, 80.0,     -1,    -1 },
	{ Geo::getLayerID( 3, 2),    1.,  5.,  5.,   0.001,  0.4 }
      }
    };

  CrCutsPtOld11MiddleRad = 
    {
       1, // use vertex
     -0.9, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 3, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 3, 1),   0.5,  20.0,  20.0,     -1,    -1 },
	{ Geo::getLayerID( 3, 2),   0.5,  1.6,  1.6,  0.001,  0.1 }
      }
    };
 


  CrCutsPt01RadMv1 = CrCutsPt01Rad;
  {
    CrCutsPt01RadMv1.layerCuts[0].iLayer = Geo::getLayerID( 0, 1);
    CrCutsPt01RadMv1.layerCuts[0].iLayer = Geo::getLayerID( 0, 2);
   CrCutsPt01RadMv1.layerCuts[0].iLayer = Geo::getLayerID( 0, 3);
  }


 

 //==========  Forward cuts

 
  CrCutsPt1Fwd = 
    {
       1, // use vertex
     .45, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 2, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 1), 0.007, 0.09, 0.05,     -1,    -1 },
	{ Geo::getLayerID( 2, 2), 0.008, 0.05, 0.009, 0.001, 0.005 }
      }
    };


  CrCutsPt05Fwd = 
    {
       1, // use vertex
     .3, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 2, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 1),  0.02,  0.1,  0.06,     -1,    -1 },
	{ Geo::getLayerID( 2, 2),  0.02,  0.06, 0.015, 0.001, 0.006 }
      }
    };

  CrCutsPt05FwdL = 
    {
       1, // use vertex
     .45, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 2, 4),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 5),  0.02,  0.1,  0.06,     -1,    -1 },
	{ Geo::getLayerID( 2, 6),  0.02,  0.06, 0.015, 0.001, 0.006 }
      }
    };
 
 CrCutsPt05FwdC = 
    {
       1, // use vertex
    .3, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 2, 2),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 3),  0.22,  0.35,  0.35,     -1,    -1 },
	{ Geo::getLayerID( 2, 4),  0.12,  0.06, 0.16, 0.001, 0.006 }
      }
    };

  CrCutsPt05FwdMiddle = 
    {
       1, // use vertex
    .45, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 5, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 5, 1),  0.06,  1.0,  0.5,     -1,    -1 },
	{ Geo::getLayerID( 5, 2),  0.06,  0.4, 0.22,  0.001,  0.08 }
      }
    };

  CrCutsPt05FwdMiddleL = 
    {
        1, // use vertex
    .45, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 5, 3),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 5, 4),  0.06,  3.5,  3.5,     -1,    -1 },
	{ Geo::getLayerID( 5, 5),  0.06,  0.45, 0.35,  0.001,  0.25 }
      }
    };

  CrCutsPt05FwdMiddleI1 = 
    {
       1, // use vertex
     .45, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 3, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 5, 0),   0.3, 10.0,   9.,     -1,    -1 },
	{ Geo::getLayerID( 5, 1),   0.2,  0.8, 0.6,    0.07,  0.1 }
      }
    };

  


  CrCutsPt03Fwd = 
    {
        1, // use vertex
    -.1, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 2, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 1),  0.03,  0.1,  0.06,     -1,    -1 },
	{ Geo::getLayerID( 2, 2),  0.03,  0.06, 0.25, 0.001, 0.006 }
      }
    };

  CrCutsPt01Fwd = 
    {
       1, // use vertex
     -.1, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 2, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 2, 1),  0.16,  1.2,  0.2,     -1,    -1 },
	{ Geo::getLayerID( 2, 2),  0.15,  0.5,  0.25, 0.001, 0.02 }
      }
    };



  CrCutsPt1MiddleFwd = 
    {
       1, // use vertex
     0.5, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 5, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 5, 1), 0.007, 0.09, 0.05,     -1,    -1 },
	{ Geo::getLayerID( 5, 2), 0.008, 0.05, 0.009, 0.001, 0.005 }
      }
    };


  CrCutsPt01MiddleFwd = 
    {
       1, // use vertex
     -0.5, // pt
      { //    layer           phiL    TL    TH      dv    duz
	{ Geo::getLayerID( 5, 0),    -1,   -1,   -1,     -1,    -1 },
	{ Geo::getLayerID( 5, 1), 0.2, 1.2, 0.2,     -1,    -1 },
	{ Geo::getLayerID( 5, 2), 0.2, 0.5, 0.25, 0.001, 0.03 }
      }
    }; 





  // ===========   prolongation cuts ====

  //  cuts for fwd/bckwd pt>=1
  PrCutsPt1Fwd = {
    {  // phi       T        V       UZ      dupV      dupUZ
      { 10*0.002, 5*0.4, 1.5*0.05, 10*0.03, 1*0.0050, 5*0.0050 }, // vol 0
      { 1*0.003, 1*0.06, 1.*0.04,  .3*0.04, 1*0.0050, 1*0.0050 }, // vol 1,2
      { 1*0.003, 1*0.06, 1.*0.04,  .3*0.04, 1*0.0050, 1*0.0050 },  
      { 2*0.005, 1*2.0,  2.*0.2,   2.*0.3,  5*0.02,   5*0.12   }, // vol 3
      { 2*0.01,  1*0.5,  1.*0.5,   .5*0.5,  1*0.02,   5*0.02   }, // vol 4,5
      { 2*0.01,  1*0.5,  1.*0.5,   .5*0.5,  1*0.02,   5*0.02   }, 
      { 1*0.005, 1*2.0,  2*0.3,    2.*1.5,  5*0.1,    5*1.1    }, // vol 6
      { 1*0.02,  2*1.0,  2*2.0,    2.*1.2,  5*0.2,    5*0.2    }, // vol 7,8
      { 1*0.02,  2*1.0,  2*2.0,    2.*1.2,  5*0.2,    5*0.2    }
    }
  };

  //  cuts for fwd/bckwd pt>=0.5
  PrCutsPt05Fwd = {
    {  // phi    T       V     UZ     dupV   dupUZ
      { 0.02,    2.,   0.08,  0.32,  0.007,  0.035  }, // vol 0
      { 0.004, 0.06,  0.042,  0.03,  0.005,  0.004  }, // vol 1,2
      { 0.004, 0.06,  0.042,  0.03,  0.005,  0.004  },  
      { 0.02,   2.0,    0.4,   0.3,   0.02,   0.12  }, // vol 3 - empty
      { 0.035,  0.5,    0.9,  0.61,  0.015,    0.1  }, // vol 4,5
      { 0.035,  0.5,    0.9,  0.61,  0.015,    0.1  }, 
      { 0.005,  2.0,    0.6,    3.,    0.5,    5.0  }, // vol 6 - empty
      { 0.02,   2.0,    4.0,    3.,     1.,    1.0  }, // vol 7,8 - empty
      { 0.02,   2.0,    4.0,    3.,     1.,    1.0  }
    }
  };

  PrCutsPt05FwdMiddle = {
    {  // phi    T       V     UZ     dupV   dupUZ
      { 0.03,  0.8,   0.06,   0.5,  0.005,  0.025  }, // vol 0
      { 0.02,  0.22,   0.3,  0.35,  0.005,  0.002  }, // vol 1,2
      { 0.02,  0.22,   0.3,  0.35,  0.005,  0.002  },  
      { 0.01,   2.0,  0.15,  0.35,  0.005,  0.035  }, // vol 3
      { 0.008,  0.41,  0.4,   0.4,   0.02,    0.1  }, // vol 4,5
      { 0.008,  0.41,  0.4,   0.4,   0.02,    0.1  }, 
      { 0.005,  2.0,   0.6,    3.,    0.5,    5.0  }, // vol 6 - empty
      { 0.03,   1.5,   2.0,   1.8,    0.3,    0.8  }, // vol 7,8 
      { 0.03,   1.5,   2.0,   1.8,    0.3,    0.8  }
    }
  };

  PrCutsPt05FwdMiddleL = {
    {  // phi    T       V     UZ     dupV   dupUZ
      { 0.03,  1.6,   0.06,   0.4,  0.007,  0.035  }, // vol 0
      { 0.03,  0.5,   0.5,  0.4,  0.004,  0.002  }, // vol 1,2
      { 0.03,  0.5,   0.5,  0.4,  0.004,  0.002  },  
      { 0.01,   2.0,  0.15,  0.35,  0.005,  0.035  }, // vol 3 - empty
      { 0.015,  0.7,  0.55,   0.5,   0.02,    0.08  }, // vol 4,5
      { 0.015,  0.7,  0.55,   0.5,   0.02,    0.08  }, 
      { 0.005,  2.0,   0.6,    3.,    0.5,    5.0  }, // vol 6 - empty
      { 0.03,   1.5,   2.0,   1.8,    0.3,    0.8  }, // vol 7,8 - empty
      { 0.03,   1.5,   2.0,   1.8,    0.3,    0.8  }
    }
  };

  PrCutsPt05FwdMiddleI1 = {
    {  // phi    T       V     UZ     dupV   dupUZ
      { 0.08,  0.1,   0.12,   0.3,  0.21,  0.025  }, // vol 0
      { 0.02,  0.22,   0.1,  0.1,  0.01,  0.002  }, // vol 1,2
      { 0.02,  0.22,   0.1,  0.1,  0.01,  0.002  },  
      { 0.012, 2.0,   0.25,  0.1,  0.055,  0.035  }, // vol 3
      { 0.008,  0.41,  0.4,   0.2,   0.08,    0.05  }, // vol 4,5
      { 0.008,  0.41,  0.4,   0.2,   0.08,    0.05  }, 
      { 0.005,  2.0,   0.6,    3.,    0.5,    5.0  }, // vol 6 - empty
      { 0.03,   1.5,   2.0,   1.8,    0.3,    0.8  }, // vol 7,8 
      { 0.03,   1.5,   2.0,   1.8,    0.3,    0.8  }
    }
  };


  PrCutsPt01Fwd = {
    {  // phi    T       V     UZ     dupV   dupUZ
      { 0.11,  2.1,    0.4,  0.8,  0.007,  0.04  }, // vol 0
      { 0.02, 0.25,    0.2,  0.2,  0.005,  0.004  }, // vol 1,2
      { 0.02, 0.25,    0.2,  0.2,  0.005,  0.004  },  
      { 0.02,   2.0,    0.4,  0.3,   0.02,   0.12  }, // vol 3 - empty
      { 0.08,  1.5,    2.0,   2.5,  0.02,    0.1  }, // vol 4,5
      { 0.08,  1.5,    2.0,   2.5,  0.02,    0.1  }, 
      { 0.005,  2.0,    0.6,    3.,    0.5,    5.0  }, // vol 6 - empty
      { 0.02,   2.0,    4.0,    3.,     1.,    1.0  }, // vol 7,8 - empty
      { 0.02,   2.0,    4.0,    3.,     1.,    1.0  }
    }
  };

 

  //  cuts for radial pt>=0.5 very loose cuts
 PrCutsPt05RadLoose = { 
  {  // phi  T       V     UZ    dupV    dupUZ    
    { 0.005, 0.5,  0.06, 0.15, 5*0.0050, 0.025 }, // vol 0
    { 0.006, 0.2,  0.05, 0.04, 5*0.0050, 0.005 }, // vol 1,2
    { 0.006, 0.2,  0.05, 0.04, 5*0.0050, 0.005 },  
    { 0.02,  2.5,  0.5,  0.5,  5*0.02,   0.11   }, // vol 3
    { 0.03,  0.5,  0.9,  0.7,  5*0.02,   0.12   }, // vol 4,5
    { 0.03,  0.5,  0.9,  0.7,  5*0.02,   0.12   }, 
    { 0.01,  1.6,  0.5,  1.5,  5*0.1,    1.1    }, // vol 6
    { 0.03,  1.6,  1.1,  1.6,  5*0.2,    0.9    }, // vol 7,8
    { 0.03,  1.6,  1.1,  1.6,  5*0.2,    0.9    }
  }
};

//  cuts for radial pt>=0.5 
 PrCutsPt05Rad = { 
  {  // phi  T       V     UZ    dupV    dupUZ    
    { 0.005, 0.5,  0.06, 0.05, 5*0.0050, 0.025 }, // vol 0
    { 0.006, 0.2,  0.05, 0.02, 5*0.0050, 0.005 }, // vol 1,2
    { 0.006, 0.2,  0.03, 0.02, 5*0.0050, 0.005 },  
    { 0.02,  2.5,  0.3,  0.35,  5*0.02,   0.11   }, // vol 3
    { 0.03,  0.5,  0.4,  0.3,  5*0.02,   0.12   }, // vol 4,5
    { 0.03,  0.5,  0.4,  0.3,  5*0.02,   0.12   }, 
    { 0.01,  1.6,  0.3,  1.5,  5*0.1,    1.1    }, // vol 6
    { 0.03,  1.6,   2.,  1.5,  5*0.2,    0.9    }, // vol 7,8
    { 0.03,  1.6,   2.,  1.5,  5*0.2,    0.9    }
  }
};

 PrCutsPt05RadMiddle = { 
  {  // phi  T       V     UZ    dupV    dupUZ    
    { 0.006, 0.31, 0.06, 0.15,  0.005, 0.004 }, // vol 0
    { 0.006, 0.2,  0.06, 0.02,  0.005, 0.005 }, // vol 1,2 -empty
    { 0.006, 0.2,  0.06, 0.02,  0.005, 0.005 },  
    { 0.006, 1.0,  0.2,  0.25,  0.008,  0.1   }, // vol 3
    { 0.02,  0.3,  0.07,  0.07, 0.025,  0.02   }, // vol 4,5
    { 0.02,  0.3,  0.07,  0.07, 0.025,  0.02   }, 
    { 0.01,  1.5,  0.21,  1.2,  0.11,   1.0    }, // vol 6
    { 0.01,  1.2,    1.,  1.1,  0.2,    0.8    }, // vol 7,8
    { 0.01,  1.2,    1.,  1.1,  0.2,    0.8    }
  }
};


 
 

 //  cuts for radial pt>=0.3
 PrCutsPt03Rad = { 
   {  // phi  T       V      UZ    dupV    dupUZ    
     { 0.012, 0.5,  0.1,     0.2,  0.0050, 10*0.005 }, // vol 0 
     { 0.01,  0.13,  0.085, 0.09,  0.0050,  0.0015 }, // vol 1,2
     { 0.01,  0.13,  0.085, 0.09,  0.0050,  0.0015 },  
     { 0.02,  2.6,  0.6,     0.6,   0.015,  0.2   }, // vol 3
     { 0.015,  0.8,  0.8,    0.6,   0.03,     0.1   }, // vol 4,5
     { 0.015,  0.8,  0.8,    0.6,   0.03,     0.1   }, 
     { 0.03,   4.0,  2.,     1.2,   0.04,   1.0    }, // vol 6
     { 0.05,   2.5,  4.,     2.0,   0.5,    3.5    }, // vol 7,8
     { 0.05,   2.5,  4.,     2.0,   0.5,    3.5    }
   }
 }; 

 //  cuts for radial pt>=0.1 
 PrCutsPt01Rad = { 
   {  // phi    T     V     UZ    dupV    dupUZ    
     { 0.02,   1.0,  0.4,  0.35,  0.008, 0.015 }, // vol 0
     { 0.009, 0.18,  0.06,  0.12, 0.002, 0.005 }, // vol 1,2
     { 0.009, 0.18,  0.06,  0.12, 0.002, 0.005 },  
     { 0.1,    3.0,  1.1,  1.2,   0.12,   0.2   }, // vol 3
     { 0.05,    .6,  1.4,  1.1,   0.08,   0.24   }, // vol 4,5
     { 0.05,    .6,  1.4,  1.1,   0.08,   0.24   }, 
     { 0.055,  2.0,  1.5,  1.2,   0.05,   1.0    }, // vol 6
     { 0.2,    3.0,  3.0,  2.0,    1.5,    2.0    }, // vol 7,8
     { 0.2,    3.0,  3.0,  2.0,    1.5,    2.0    }
   }
 };



  //  cuts for radial pt>=0.1 
 PrCutsPt01MiddleRad = { 
   {  // phi  T       V     UZ    dupV    dupUZ    
     { 0.02,  1.0,  0.4,  0.8,  0.005, 0.01 }, // vol 0
     { 0.15,  1.,   1.5,  1.3,  0.005, 0.005 }, // vol 1,2
     { 0.15,  1.,   1.5,  1.3,  0.005, 0.005 },  
     { 0.05,  2.0,  2.0,  2.0,  0.02,   0.1   }, // vol 3
     { 0.2,   4.0,  4.0,  5.,   0.05,   0.1   }, // vol 4,5
     { 0.2,   4.0,  4.0,  5.,   0.05,   0.1   }, 
     { 0.4,   4.0,  4.0,  3.0,   0.1,   0.5    }, // vol 6
     { 0.2,   4.0,  10.0,  4.0,  0.5,   0.5    }, // vol 7,8
     { 0.2,   4.0,  10.0,  4.0,   0.5,  0.5    }
   }
 };

}
