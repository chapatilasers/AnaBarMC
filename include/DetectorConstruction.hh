#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"

//---------------------------------------------------------------------------

class G4VPhysicalVolume;
class DetectorSD;
class DetectorMessenger;

//---------------------------------------------------------------------------

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();

  void UpdateGeometry();

  inline G4VPhysicalVolume* GetExpHall()    { return fExpHall;  };
  inline G4VPhysicalVolume* GetDet1Vol()    { return fDet1Vol;  };
  inline G4VPhysicalVolume* GetDet2Vol()    { return fDet2Vol;  };
  inline G4VPhysicalVolume* GetFinger()    { return fFinger;  };
  inline G4VPhysicalVolume* GetAnaBar()    { return fAnaBar;  };
  inline DetectorSD*        GetDetSD()      { return fDetSD;    };

  inline void SetTumourOn     ( G4int    tumon )   { fTumourOn     = tumon; }
  inline void SetTumourRadius ( G4double radius )  { fTumourRadius = radius; }
  inline void SetTumourHeight ( G4double height )  { fTumourHeight = height; }

  private:

  G4NistManager*     fNistManager;
  DetectorMessenger* fDetMessenger;

  G4VPhysicalVolume* fExpHall;
  G4VPhysicalVolume* fDet1Vol;
  G4VPhysicalVolume* fDet2Vol;
  G4VPhysicalVolume* fFinger;
  G4VPhysicalVolume* fAnaBar;

  DetectorSD*        fDetSD;

  G4int              fTumourOn;
  G4double           fTumourRadius;
  G4double           fTumourHeight;

};
#endif

//---------------------------------------------------------------------------

