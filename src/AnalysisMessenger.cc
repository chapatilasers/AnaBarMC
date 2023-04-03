#include "AnalysisMessenger.hh"
#include "AnalysisManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"

//---------------------------------------------------------------------------

AnalysisMessenger::AnalysisMessenger(AnalysisManager* anMana)
:fAnalysisManager(anMana)
{
  fAnalysisDir = new G4UIdirectory("/AnaBarMC/analysis/");
  fAnalysisDir->SetGuidance("Analysis output control");

  fOutFileCmd = new G4UIcmdWithAString("/AnaBarMC/analysis/setOutputFile",this);
  fOutFileCmd->SetGuidance("Set the full name and path of the output ROOT file");
  fOutFileCmd->SetParameterName("choice",true);
  fOutFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fGeoFileCmd = new G4UIcmdWithAString("/AnaBarMC/analysis/setGeometryFile",this);
  fGeoFileCmd->SetGuidance("Set the full name and path of the output Geometry file");
  fGeoFileCmd->SetParameterName("choice",true);
  fGeoFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
} 

//---------------------------------------------------------------------------

AnalysisMessenger::~AnalysisMessenger()
{
  delete fOutFileCmd;
  delete fGeoFileCmd;
}

//---------------------------------------------------------------------------

void AnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == fOutFileCmd)
    {fAnalysisManager->SetOutFileName(newValue.data());}
  if(command == fGeoFileCmd)
    {fAnalysisManager->SetGeoFileName(newValue.data());}
}

//---------------------------------------------------------------------------
