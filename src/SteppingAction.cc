#include "SteppingAction.hh"
#include "PMTSD.hh"

#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpBoundaryProcess.hh"

//---------------------------------------------------------------------------

SteppingAction::SteppingAction()
  :fOneStepPrimaries(false)
{}

//---------------------------------------------------------------------------

SteppingAction::~SteppingAction()
{}

//---------------------------------------------------------------------------

void SteppingAction::UserSteppingAction(const G4Step * theStep)
{

  G4Track* theTrack = theStep->GetTrack();
  
  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  static G4OpBoundaryProcess* boundary=NULL;

  //find the boundary process only once
  if(!boundary){
    G4ProcessManager* pm 
      = theStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
	boundary = (G4OpBoundaryProcess*)(*pv)[i];
	break;
      }
    }
  }

  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()) {
    
    boundaryStatus=boundary->GetStatus();
    
    switch(boundaryStatus){
    case Absorption:
      break;
    case Detection:	              
      {
	//Triger sensitive detector manually since photon is
	//absorbed but status was Detection.
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4String sdName="PMTSD";
	PMTSD* detSD = (PMTSD*)SDman
	  ->FindSensitiveDetector(sdName);
	if(detSD)
	  detSD->ProcessHits_constStep(theStep,NULL);
	break;
      }
    case FresnelReflection:
    case TotalInternalReflection:
    case SpikeReflection:
      break;
    default:
      break;
    }
  }
}

//---------------------------------------------------------------------------


















