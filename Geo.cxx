#include "Geo.h"
#include <iostream>

 constexpr int Geo::NVolumes;
 constexpr int Geo::NLayers; // 48

 Volume Geo::volumes[NVolumes];
 Layer  Geo::layers[NLayers];
 


void Volume::init( int s, int t, int nL, int &nLayersTot)
{
  size = s;
  type = t;
  nLayers = nL;
  for( int i=0; i<nLayers; i++ ) layerIDs[i] = nLayersTot++;
}


void Geo::init()
{
  // Volume:
  //  size; // radial position: 0 = inner, 1 = middle, 2 = outer 
  //  type; // 0 = radial modules, 1 = vertical modules
  //  nLayers; // number of layers [2..7]

  for( int i=0; i<NLayers; i++ ){
    Layer &l = layers[i];
    l.id = i; // general layer id
  }

  int nLayersTot = 0;
  volumes[0].init( 0,  0, 4, nLayersTot);
  volumes[1].init( 0, -1, 7, nLayersTot);
  volumes[2].init( 0,  1, 7, nLayersTot);

  volumes[3].init( 1,  0, 4, nLayersTot);
  volumes[4].init( 1, -1, 6, nLayersTot);
  volumes[5].init( 1,  1, 6, nLayersTot);

  volumes[6].init( 2,  0, 2, nLayersTot);
  volumes[7].init( 2, -1, 6, nLayersTot);
  volumes[8].init( 2,  1, 6, nLayersTot);
  
  if( nLayersTot != NLayers ){
    cout<<"wrong n layers!!"<<endl;
    exit(0);
  }
  
  for( int iv=0; iv<NVolumes; iv++ ){
    Volume &v = volumes[iv];
    for( int il=0; il<v.nLayers; il++ ){      
      Layer &l = layers[v.layerIDs[il]];      
      l.volume = iv;
      l.layer = il;
      l.size = v.size;
      l.type = v.type;
    }
  }

  // read layer sizes
  {
    ifstream in("geoLayerSizes.txt");    
    if( !in.is_open() ){
      cout<<"Geo:: Can not open geoLayerSizes.txt file"<<endl;
      exit(0);
    }
    for( int i=0; i<NLayers; i++){
      Layer &l = layers[i];
      int j;
      in>>j>>l.r>>l.z>>l.tMin>>l.tMax;
      if( j!=i ){
	cout<<"Geo: geoLayerSizes.txt file broken"<<endl;
	exit(1);
      }
    }
    in.close();
  }

   // read layer field
  {
   static constexpr double CLight = 0.000299792458; // speed of light 
   ifstream in("geoLayerField.txt");    
    if( !in.is_open() ){
      cout<<"Geo:: Can not open geoLayerField.txt file"<<endl;
      exit(0);
    }

    for( int il=0; il<NLayers; il++){
      Layer &l = layers[il];
      double sigma[3]={0,0,0};
      for( int ir=0; ir<3; ir++){
	int jl;
	in>>jl;
	if( jl!=il ){
	  cout<<"Geo:: geo field file broken"<<endl;
	  exit(1);
	}
	for( int i=0; i<Layer::NFieldPar; i++ ){
	  in>>l.field[ir][i];
	  l.field[ir][i]*=CLight;
	}
	in>>sigma[ir];
      }
      for( int ir=0; ir<2; ir++){
	if( sigma[ir]<=1.e-4 ){
	  for( int i=0; i<Layer::NFieldPar; i++ ) l.field[ir][i]=l.field[2][i];	
	}
      }
    }
    in.close();
  }


  for( int i=0; i<NLayers; i++ ){
    //Layer &l = layers[i];
    //cout<<"layer id "<<i<<" vol "<<l.volume<<" il "<<l.layer<<endl;
  }

}
