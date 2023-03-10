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
  fAnaBarYposCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/AnaBarYpos",this);
  fAnaBarYposCmd->SetGuidance("Set Anabar x position in cm.");
  fAnaBarZposCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/AnaBarZpos",this);
  fAnaBarZposCmd->SetGuidance("Set Anabar x position in cm.");
  fAnaBarLengthCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/AnaBarLength",this);
  fAnaBarLengthCmd->SetGuidance("Set Anabar Length in cm.");
  
  fNumberOfLayersCmd = new G4UIcmdWithAnInteger("/AnaBarMC/detector/NumberOfLayers",this);
  fNumberOfLayersCmd->SetGuidance("Set Number of AnaBar Layers.");
  fNumberOfBarsCmd = new G4UIcmdWithAnInteger("/AnaBarMC/detector/NumberOfBars",this);
  fNumberOfBarsCmd->SetGuidance("Set Number of AnaBar Bars.");
  fNumberOfSidesCmd = new G4UIcmdWithAnInteger("/AnaBarMC/detector/NumberOfSides",this);
  fNumberOfSidesCmd->SetGuidance("Set Number of AnaBar Sides.");
  fNumberOfModulesCmd = new G4UIcmdWithAnInteger("/AnaBarMC/detector/NumberOfModules",this);
  fNumberOfModulesCmd->SetGuidance("Set Number of AnaBar Modules.");
  fNumberOfPlanesCmd = new G4UIcmdWithAnInteger("/AnaBarMC/detector/NumberOfPlanes",this);
  fNumberOfPlanesCmd->SetGuidance("Set Number of AnaBar Planes.");
  
  fAnaBarWidthCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/AnaBarWidth",this);
  fAnaBarWidthCmd->SetGuidance("Set Anabar Width in cm.");
  fAnaBarThicknessCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/AnaBarThickness",this);
  fAnaBarThicknessCmd->SetGuidance("Set Anabar Thickness in cm.");
  fFingerLengthCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/FingerLength",this);
  fFingerLengthCmd->SetGuidance("Set Finger Length in cm.");
  fFingerWidthCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/FingerWidth",this);
  fFingerWidthCmd->SetGuidance("Set Finger Width in cm.");
  fFingerThicknessCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/FingerThickness",this);
  fFingerThicknessCmd->SetGuidance("Set Finger Thickness in cm.");
  fFibreDiameterCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/FibreDiameter",this);
  fFibreDiameterCmd->SetGuidance("Set Fibre Diameter in cm.");
  fFibreLengthCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/FibreLength",this);
  fFibreLengthCmd->SetGuidance("Set Fibre Length in cm.");
  fHoleDiameterCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/HoleDiameter",this);
  fHoleDiameterCmd->SetGuidance("Set Hole Diameter in cm.");
  fHoleLengthCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/HoleLength",this);
  fHoleLengthCmd->SetGuidance("Set Hole Length in cm.");
  fCladdingDiameterCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/CladdingDiameter",this);
  fCladdingDiameterCmd->SetGuidance("Set Cladding Diameter in cm.");
  fCladdingLengthCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/CladdingLength",this);
  fCladdingLengthCmd->SetGuidance("Set Cladding Length in cm.");
  fPhotoCathodeDiameterCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/PhotoCathodeDiameter",this);
  fPhotoCathodeDiameterCmd->SetGuidance("Set PhotoCathode Diameter in cm.");
  fPhotoCathodeThicknessCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/PhotoCathodeThickness",this);
  fPhotoCathodeThicknessCmd->SetGuidance("Set PhotoCathode Thickness in cm.");
  fMirrorThicknessCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/MirrorThickness",this);
  fMirrorThicknessCmd->SetGuidance("Set Mirror Thickness in cm.");
  fMylarThicknessCmd = new G4UIcmdWithADouble("/AnaBarMC/detector/MylarThickness",this);
  fMylarThicknessCmd->SetGuidance("Set Mylar Thickness in cm.");

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
  delete fAnaBarYposCmd;
  delete fAnaBarZposCmd;
  delete fNumberOfLayersCmd;
  delete fNumberOfBarsCmd;
  delete fNumberOfSidesCmd;
  delete fNumberOfModulesCmd;
  delete fNumberOfPlanesCmd;
  delete fAnaBarLengthCmd;
  delete fAnaBarWidthCmd;
  delete fAnaBarThicknessCmd;
  delete fFingerLengthCmd;
  delete fFingerWidthCmd;
  delete fFingerThicknessCmd;
  delete fFibreDiameterCmd;
  delete fFibreLengthCmd;
  delete fHoleDiameterCmd;
  delete fHoleLengthCmd;
  delete fCladdingDiameterCmd;
  delete fCladdingLengthCmd;
  delete fPhotoCathodeDiameterCmd;
  delete fPhotoCathodeThicknessCmd;
  delete fMirrorThicknessCmd;
  delete fMylarThicknessCmd;

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
  if (command == fAnaBarYposCmd)
    fDetector->SetAnaBarYpos(fAnaBarYposCmd->GetNewDoubleValue(newValue));
  if (command == fAnaBarZposCmd)
    fDetector->SetAnaBarZpos(fAnaBarZposCmd->GetNewDoubleValue(newValue));
  
  if (command == fNumberOfLayersCmd)
    fDetector->SetNumberOfLayers(fNumberOfLayersCmd->GetNewIntValue(newValue));
  if (command == fNumberOfBarsCmd)
    fDetector->SetNumberOfBars(fNumberOfBarsCmd->GetNewIntValue(newValue));
  if (command == fNumberOfSidesCmd)
    fDetector->SetNumberOfSides(fNumberOfSidesCmd->GetNewIntValue(newValue));
  if (command == fNumberOfModulesCmd)
    fDetector->SetNumberOfModules(fNumberOfModulesCmd->GetNewIntValue(newValue));
  if (command == fNumberOfPlanesCmd)
    fDetector->SetNumberOfPlanes(fNumberOfPlanesCmd->GetNewIntValue(newValue));
  
  if (command == fAnaBarLengthCmd)
    fDetector->SetAnaBarLength(fAnaBarLengthCmd->GetNewDoubleValue(newValue));
  if (command == fAnaBarWidthCmd)
    fDetector->SetAnaBarWidth(fAnaBarWidthCmd->GetNewDoubleValue(newValue));
  if (command == fAnaBarThicknessCmd)
    fDetector->SetAnaBarThickness(fAnaBarThicknessCmd->GetNewDoubleValue(newValue));
  
  if (command == fFingerLengthCmd)
    fDetector->SetFingerLength(fFingerLengthCmd->GetNewDoubleValue(newValue));
  if (command == fFingerWidthCmd)
    fDetector->SetFingerWidth(fFingerWidthCmd->GetNewDoubleValue(newValue));
  if (command == fFingerThicknessCmd)
    fDetector->SetFingerThickness(fFingerThicknessCmd->GetNewDoubleValue(newValue));
  
  if (command == fFibreDiameterCmd)
    fDetector->SetFibreDiameter(fFibreDiameterCmd->GetNewDoubleValue(newValue));
  if (command == fFibreLengthCmd)
    fDetector->SetFibreLength(fFibreLengthCmd->GetNewDoubleValue(newValue));
  
  if (command == fHoleDiameterCmd)
    fDetector->SetHoleDiameter(fHoleDiameterCmd->GetNewDoubleValue(newValue));
  if (command == fHoleLengthCmd)
    fDetector->SetHoleLength(fHoleLengthCmd->GetNewDoubleValue(newValue));
  
  if (command == fCladdingDiameterCmd)
    fDetector->SetCladdingDiameter(fCladdingDiameterCmd->GetNewDoubleValue(newValue));
  if (command == fCladdingLengthCmd)
    fDetector->SetCladdingLength(fCladdingLengthCmd->GetNewDoubleValue(newValue));
  
  if (command == fPhotoCathodeThicknessCmd)
    fDetector->SetPhotoCathodeThickness(fPhotoCathodeThicknessCmd->GetNewDoubleValue(newValue));
  if (command == fPhotoCathodeDiameterCmd)
    fDetector->SetPhotoCathodeDiameter(fPhotoCathodeDiameterCmd->GetNewDoubleValue(newValue));
  if (command == fMirrorThicknessCmd)
    fDetector->SetMirrorThickness(fMirrorThicknessCmd->GetNewDoubleValue(newValue));
  if (command == fMylarThicknessCmd)
    fDetector->SetMylarThickness(fMylarThicknessCmd->GetNewDoubleValue(newValue));
}

//---------------------------------------------------------------------------
