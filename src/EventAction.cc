#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"

#include "EventAction.hh"
#include "DetectorHit.hh"
#include "PMTHit.hh"
#include "AnalysisManager.hh"
#include "PrimaryGeneratorAction.hh"

//---------------------------------------------------------------------------

EventAction::EventAction( AnalysisManager* ana, PrimaryGeneratorAction* pga )
  :fAnaManager(ana), fPGA(pga), fDetCollID(-1), fPMTCollID(-1)
{;}

//---------------------------------------------------------------------------

EventAction::~EventAction()
{;}

//---------------------------------------------------------------------------

void EventAction::BeginOfEventAction(const G4Event* evt)
{ 
  if( evt->GetEventID() == 0 )
    fAnaManager->InitOutput();
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(fDetCollID < 0 ) 
    fDetCollID = SDman->GetCollectionID("DetSDCollection");
  if(fPMTCollID < 0 ) 
    fPMTCollID = SDman->GetCollectionID("PMTSDCollection");

}

//---------------------------------------------------------------------------

void EventAction::EndOfEventAction(const G4Event* evt)
{

  G4int event_id   = evt->GetEventID();  
  if ( event_id%1 == 0 ) 
    G4cout <<"Event " << event_id << G4endl;

  fAnaManager->ZeroArray(); 

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  DetectorHitsCollection *DHC = NULL;
  PMTHitsCollection      *PHC = NULL;
  G4int det_hits = 0;
  G4int pmt_hits = 0;

  if(HCE) {
    if( fDetCollID >= 0 ) DHC = (DetectorHitsCollection*)(HCE->GetHC(fDetCollID));       
    if( fPMTCollID >= 0 ) PHC = (PMTHitsCollection*)(HCE->GetHC(fPMTCollID));
  }

  // Detector
  if(DHC) {
    det_hits = DHC->entries();
    if( det_hits != 0 ) {
      for( G4int i = 0; i < det_hits; i++) {
	DetectorHit* hit1 = static_cast<DetectorHit*>( DHC->GetHit(i) );
	fAnaManager->SetStepPDef( (G4ParticleDefinition*) hit1->GetPDef() );
	fAnaManager->SetStepP3( (G4ThreeVector) hit1->GetMom() );
	fAnaManager->SetStepPosPre( (G4ThreeVector) hit1->GetPosPre() );
	fAnaManager->SetStepPosPost( (G4ThreeVector) hit1->GetPosPost() );
	fAnaManager->SetStepTime( (G4double) hit1->GetTime() );
	fAnaManager->SetStepEdep( (G4double) hit1->GetEdep() );
	fAnaManager->SetStepID( (G4int) hit1->GetID() );
	fAnaManager->FillArray( i ); 
      }
    }
  }

  // PMT
  if(PHC) {
    pmt_hits = PHC->entries();
    //std::cout << "Number of PMT Hits at end of event = " << pmt_hits << std::endl;
    if( pmt_hits != 0 ) {
      for( G4int j = 0; j < pmt_hits; j++) {
	PMTHit* hit2 = static_cast<PMTHit*>( PHC->GetHit(j) );
	if (hit2->GetPMTNumber() == 14) {
		fAnaManager->SetPhotonCountZero( (G4int) hit2->GetPhotonCount() );
	}else{
		if (hit2->GetPMTNumber() == 0) {
			fAnaManager->SetPhotonCountOne( (G4int) hit2->GetPhotonCount() );
		}
		if (hit2->GetPMTNumber() == 1) {
			fAnaManager->SetPhotonCountTwo( (G4int) hit2->GetPhotonCount() );
		}
		if (hit2->GetPMTNumber() == 2) {
			fAnaManager->SetPhotonCountThree( (G4int) hit2->GetPhotonCount() );
		}
	}
	fAnaManager->SetPMTNumber( (G4int) hit2->GetPMTNumber() );
	//std::cout << "hit " << j << " pmt number " << hit2->GetPMTNumber() << " number of photons = " << hit2->GetPhotonCount() << std::endl;
      }
    }
  }
  
  // Primary
  if( det_hits != 0 ) {
    fAnaManager->SetPrimaryDirection ( (G4ThreeVector)fPGA->GetDirection() );
    fAnaManager->SetPrimaryEnergy    ( (G4double)fPGA->GetEnergy() );
    fAnaManager->SetPrimaryTime      ( (G4double)fPGA->GetTime() );
    fAnaManager->SetPrimaryPDef      ( (G4ParticleDefinition*)fPGA->GetPrimPDef() );
    
    fAnaManager->FillTree(); 
  }

}

//---------------------------------------------------------------------------
