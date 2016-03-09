#include "PMTHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

//--------------------------------------------------------------------------

G4Allocator<PMTHit> PMTHitAllocator;

PMTHit::PMTHit()
{
  pmtNumber = -1;
  photons   = 0;
}

//--------------------------------------------------------------------------

PMTHit::~PMTHit()
{;}

//--------------------------------------------------------------------------

PMTHit::PMTHit(const PMTHit &right)
  : G4VHit()
{
  pmtNumber = right.pmtNumber;
  photons   = right.photons;
}

//--------------------------------------------------------------------------

const PMTHit& PMTHit::operator=(const PMTHit &right)
{
  pmtNumber = right.pmtNumber;
  photons   = right.photons;
  return *this;
}

//--------------------------------------------------------------------------

G4int PMTHit::operator==(const PMTHit&) const
{ return 0;}

//--------------------------------------------------------------------------

void PMTHit::Draw()
{;}

//--------------------------------------------------------------------------

void PMTHit::Print()
{;}

//--------------------------------------------------------------------------








