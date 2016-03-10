#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

//---------------------------------------------------------------------------

DetectorMessenger::DetectorMessenger(DetectorConstruction* Detect)
  :fDetector(Detect)
{
  fDetectorDir        = new G4UIdirectory("/AnaBarMC/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");
  
  fTumourOnCmd        = new G4UIcmdWithAnInteger("/AnaBarMC/detector/TumourOn",this);
  fTumourOnCmd->SetGuidance("Set 1 for tumour on, 0 for off.");

  fTumourRadiusCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/TumourRadius",this);
  fTumourRadiusCmd->SetGuidance("Set the radius of the tumour in cm.");

  fTumourHeightCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/TumourHeight",this);
  fTumourHeightCmd->SetGuidance("Set the radius of the tumour in cm.");

  fAnaBarXposCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/AnaBarXpos",this);
  fAnaBarXposCmd->SetGuidance("Set Anabar x position in cm.");

  fUpdateCmd          = new G4UIcommand("/AnaBarMC/detector/update",this);
  fUpdateCmd->SetGuidance("Update the detector geometry with changed values.");
  fUpdateCmd->SetGuidance("Must be run before beamOn if detector has been changed.");  
}

//---------------------------------------------------------------------------

DetectorMessenger::~DetectorMessenger()
{
  delete fDetectorDir;

  delete fTumourOnCmd;
  delete fTumourRadiusCmd;
  delete fTumourHeightCmd;
  delete fAnaBarXposCmd;

  delete fUpdateCmd;
}

//---------------------------------------------------------------------------

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if(command == fUpdateCmd)
    fDetector->UpdateGeometry();
  if (command == fTumourOnCmd)
    fDetector->SetTumourOn(fTumourOnCmd->GetNewIntValue(newValue));
  if (command == fTumourRadiusCmd)
    fDetector->SetTumourRadius(fTumourRadiusCmd->GetNewDoubleValue(newValue));
  if (command == fTumourHeightCmd)
    fDetector->SetTumourHeight(fTumourHeightCmd->GetNewDoubleValue(newValue));
  if (command == fAnaBarXposCmd)
    fDetector->SetAnaBarXpos(fAnaBarXposCmd->GetNewDoubleValue(newValue));
}

//---------------------------------------------------------------------------
