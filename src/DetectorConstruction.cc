#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "DetectorSD.hh"

#include "G4Material.hh"
#include "G4BooleanSolid.hh"
#include "G4CSGSolid.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4UserLimits.hh"

#include "G4TransportationManager.hh"
#include "G4SDManager.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4VisAttributes.hh"
#include "G4String.hh"
#include "globals.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

DetectorConstruction::DetectorConstruction()
{
  fDetMessenger = new DetectorMessenger(this);

  fNistManager  = G4NistManager::Instance();
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4String command = "/control/execute macros/DetectorSetup.mac";
  UI->ApplyCommand(command);

}

//---------------------------------------------------------------------------

DetectorConstruction::~DetectorConstruction() 
{
}

//---------------------------------------------------------------------------

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //vacuum
  G4double vacuumDensity  = 1.e-25 *g/cm3;
  G4double pressure       = 3.e-18*pascal;
  G4double temperature    = 2.73*kelvin;
  G4Material* vacuum      = new G4Material("Galactic", 1., 1.01*g/mole,
					   vacuumDensity,kStateGas,temperature,pressure);

  //---------------------------------------------------------------------------
  // Create Experimental Hall
  //---------------------------------------------------------------------------

  G4Material* expHall_mat = fNistManager->FindOrBuildMaterial("G4_AIR");

  G4Box* expHall_box           = new G4Box("expHall_box",
					   1.5 *m, 1.5 *m, 1.5 *m );
  
  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,
						     expHall_mat,
						     "expHall_log", 0, 0, 0);
  
  fExpHall                     = new G4PVPlacement(0, G4ThreeVector(),
						   expHall_log, "expHall", 0, false, 0);

  //---------------------------------------------------------------------------
  // Create Finger and Analyzer Bar Counters
  //---------------------------------------------------------------------------

  G4double anaBar_x = 20.0 *cm;
  G4double anaBar_y = 4.0 *cm;
  G4double anaBar_z = 4.0 *cm;
  G4ThreeVector anaBar_pos = G4ThreeVector(0.0,0.0,0.0);
  G4Material* anaBar_mat = fNistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  
  G4double finger_x = 2.54 *cm;
  G4double finger_y = 2.54 *cm;
  G4double finger_z = 0.75 *cm;
  G4ThreeVector finger_pos = G4ThreeVector(0.5*anaBar_x-0.5*finger_x,0.0,0.5*anaBar_z+0.5*finger_z);
  G4Material* finger_mat = fNistManager->FindOrBuildMaterial("G4_POLYETHYLENE");

  G4Box* finger_box           = new G4Box("finger_box",
					   finger_x/2.0, finger_y/2.0, finger_z/2.0);
  
  G4LogicalVolume* finger_log = new G4LogicalVolume(finger_box,
						     finger_mat,
						     "finger_log", 0, 0, 0);
  
  fFinger                     = new G4PVPlacement(0, finger_pos, finger_log, "finger", expHall_log, false, 0);

  G4Box* anaBar_box           = new G4Box("anaBar_box", anaBar_x/2.0, anaBar_y/2.0, anaBar_z/2.0);
  
  G4LogicalVolume* anaBar_log = new G4LogicalVolume(anaBar_box,
						     anaBar_mat,
						     "anaBar_log", 0, 0, 0);
  
  fAnaBar                     = new G4PVPlacement(0, anaBar_pos,anaBar_log, "anaBar", expHall_log, false, 1);
  
  //---------------------------------------------------------------------------
  // Set Step Limits, Sensitive Detector and Visualisation
  //---------------------------------------------------------------------------

  G4double maxStep = 0.5 *mm;;
  G4UserLimits* stepLimit = new G4UserLimits(maxStep);
  //  det1_log->SetUserLimits(stepLimit);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  fDetSD = new DetectorSD("DetSD", 2);
  SDman->AddNewDetector( fDetSD );
  anaBar_log->SetSensitiveDetector( fDetSD );
  finger_log->SetSensitiveDetector( fDetSD );
  
  G4VisAttributes* blue    = new G4VisAttributes( G4Colour(0.0,0.0,1.0)   );
  G4VisAttributes* red     = new G4VisAttributes( G4Colour(1.0,0.0,0.0)   );
  G4VisAttributes* green   = new G4VisAttributes( G4Colour(0.0,1.0,0.0)   );

  expHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  anaBar_log->SetVisAttributes(green);
  finger_log->SetVisAttributes(blue);

  return fExpHall;
}

//---------------------------------------------------------------------------

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//---------------------------------------------------------------------------
