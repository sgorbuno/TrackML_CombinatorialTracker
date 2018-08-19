#include "DataStructures.h"

using namespace std;

void Hit::Print() const
{
  cout<<"hit: vol "<<volume<<" l "<<layer<<" x "<<x<<" y "<<y<<" z "<<z<<" r "<<r<<" phi "<<phi<<endl; 
}

void HitMC::Print() const
{  
  cout<<"hitmc: "<<partID<< "  x "<<x<<" y "<<y<<" z "<<z<<" p "<<p<<" pt "<<pt<<" pz "<<pz<<" p "<<p<<endl; 
}

void  Particle::Print() const
{
  cout<<"particle: x "<<x<<" y "<<y<<" z "<<z<<" pt "<<pt<<" pz "<<pz<<" p "<<p<<endl;
};
