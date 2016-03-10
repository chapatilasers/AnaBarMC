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

  fTumourRadius = 0.5;
  fTumourHeight = 0.0;
  fAnaBarXpos   = 0.0;
  
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

  G4Box* expHall_box           = new G4Box("expHall_box",
					   1.5 *m, 1.5 *m, 1.5 *m );
  
  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,
						     fNistManager->FindOrBuildMaterial("G4_AIR"),
						     "expHall_log", 0, 0, 0);
  
  fExpHall                     = new G4PVPlacement(0, G4ThreeVector(),
						   expHall_log, "expHall", 0, false, 0);

  //---------------------------------------------------------------------------
  // Create Detectors 
  //---------------------------------------------------------------------------

   G4Material* scintillator__mat = fNistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
   G4Box* fingercounter_solid  = new G4Box("fingercounter_solid", 1.3*cm , 2.0*cm , 0.85*cm);
   G4LogicalVolume* fingercounter_log = new G4LogicalVolume(fingercounter_solid, scintillator__mat, "fingercounter_log");
   G4ThreeVector fingercounter_pos(0.0*cm , 0.0*cm , 0.0*cm);
   G4VPhysicalVolume* FingerCounter=  new G4PVPlacement(0, fingercounter_pos , fingercounter_log , "FingerCounter" , expHall_log , false , 0); 

 
   G4Box* AnaBar_solid  = new G4Box("AnaBar_solid", 11.0*cm , 2.0*cm , 2.0*cm);
   G4LogicalVolume* AnaBar_log = new G4LogicalVolume(AnaBar_solid, scintillator__mat, "AnaBar_log");
   G4ThreeVector AnaBar_pos(fAnaBarXpos*cm , 0.0*cm , -2.85*cm);
   G4VPhysicalVolume* AnaBar =  new G4PVPlacement(0, AnaBar_pos , AnaBar_log , "AnaBar" , expHall_log , false , 1); 
 


  //---------------------------------------------------------------------------
  // Set Step Limits, Sensitive Detector and Visualisation
  //---------------------------------------------------------------------------

  G4double maxStep = 0.5 *mm;;
  G4UserLimits* stepLimit = new G4UserLimits(maxStep);
  //  det1_log->SetUserLimits(stepLimit);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  fDetSD = new DetectorSD("DetSD", 2);
  SDman->AddNewDetector( fDetSD );
//  det1_log->SetSensitiveDetector( fDetSD );
    fingercounter_log->SetSensitiveDetector( fDetSD );
    AnaBar_log->SetSensitiveDetector( fDetSD );


  G4VisAttributes* blue    = new G4VisAttributes( G4Colour(0.0,0.0,1.0)   );
  G4VisAttributes* red     = new G4VisAttributes( G4Colour(1.0,0.0,0.0)   );
  G4VisAttributes* green   = new G4VisAttributes( G4Colour(0.0,1.0,0.0)   );

  expHall_log->SetVisAttributes(G4VisAttributes::Invisible);
//  det1_log->SetVisAttributes(blue);
    fingercounter_log->SetVisAttributes(blue);
    AnaBar_log->SetVisAttributes(green);
  return fExpHall;
}

//---------------------------------------------------------------------------


void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//---------------------------------------------------------------------------
