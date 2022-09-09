#ifndef PMTSD_h
#define PMTSD_h 1

#include "G4VSensitiveDetector.hh"
#include "PMTHit.hh"

class G4Step;
class G4HCofThisEvent;

class PMTSD : public G4VSensitiveDetector
{

public:
  PMTSD(G4String name, G4int Nelements);
  ~PMTSD();
  
  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  
  //A version of processHits that keeps aStep constant
  G4bool ProcessHits_constStep(const G4Step* aStep,
			       G4TouchableHistory* ROhist);
  void EndOfEvent(G4HCofThisEvent* HCE);
  
private:
  PMTHitsCollection*  fCollection;  
  G4int               fNelements;
  G4int               fHCID;
  G4int               fNhits;
  G4int*              fhitID;
  G4int*              fHits;
};

#endif

