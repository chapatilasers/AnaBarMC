#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

//---------------------------------------------------------------------------

class DetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

//---------------------------------------------------------------------------

class DetectorMessenger: public G4UImessenger
{
public:
  DetectorMessenger(DetectorConstruction*);
  ~DetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  DetectorConstruction*        fDetector;

  G4UIdirectory*               fDetectorDir;

  G4UIcmdWithAnInteger*        fTumourOnCmd;
  G4UIcmdWithADouble*          fTumourRadiusCmd;
  G4UIcmdWithADouble*          fTumourHeightCmd;
  G4UIcmdWithADouble*          fAnaBarXposCmd;
  G4UIcmdWithADouble*          fAnaBarYposCmd;
  G4UIcmdWithADouble*          fAnaBarZposCmd;
  G4UIcmdWithAnInteger*          fNumberOfLayersCmd;
  G4UIcmdWithAnInteger*          fNumberOfBarsCmd;
  G4UIcmdWithAnInteger*          fNumberOfSidesCmd;
  G4UIcmdWithAnInteger*          fNumberOfModulesCmd;
  G4UIcmdWithAnInteger*          fNumberOfPlanesCmd;
  G4UIcmdWithADouble*          fAnaBarLengthCmd;
  G4UIcmdWithADouble*          fAnaBarWidthCmd;
  G4UIcmdWithADouble*          fAnaBarThicknessCmd;
  G4UIcmdWithADouble*          fFingerLengthCmd;
  G4UIcmdWithADouble*          fFingerWidthCmd;
  G4UIcmdWithADouble*          fFingerThicknessCmd;
  G4UIcmdWithADouble*          fFingerZoffsetCmd;
  G4UIcmdWithADouble*          fFibreDiameterCmd;
  G4UIcmdWithADouble*          fFibreLengthCmd;
  G4UIcmdWithADouble*          fHoleDiameterCmd;
  G4UIcmdWithADouble*          fHoleLengthCmd;
  G4UIcmdWithADouble*          fCladdingDiameterCmd;
  G4UIcmdWithADouble*          fCladdingLengthCmd;
  G4UIcmdWithADouble*          fPhotoCathodeDiameterCmd;
  G4UIcmdWithADouble*          fPhotoCathodeThicknessCmd;
  G4UIcmdWithADouble*          fMirrorThicknessCmd;
  G4UIcmdWithADouble*          fMylarThicknessCmd;

  G4UIcommand*                 fUpdateCmd;
};

#endif

//---------------------------------------------------------------------------

