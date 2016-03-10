#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"

//---------------------------------------------------------------------------

class G4VPhysicalVolume;
class DetectorSD;
class PMTSD;
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
  inline G4VPhysicalVolume* GetFingerVol()    { return FingerCounter;  };
  inline G4VPhysicalVolume* GetAnaBarVol()    { return AnaBar;  };
  inline DetectorSD*        GetDetSD()      { return fDetSD;    };
  inline PMTSD*             GetPMTSD()      { return fPMTSD;    };

  inline void SetTumourOn     ( G4int    tumon )   { fTumourOn     = tumon; }
  inline void SetTumourRadius ( G4double radius )  { fTumourRadius = radius; }
  inline void SetTumourHeight ( G4double height )  { fTumourHeight = height; }
  inline void SetAnaBarXpos   ( G4double Xpos )    { fAnaBarXpos   = Xpos; }
 
  private:

  G4NistManager*     fNistManager;
  DetectorMessenger* fDetMessenger;

  G4VPhysicalVolume* fExpHall;
  G4VPhysicalVolume* fDet1Vol;
  G4VPhysicalVolume* fDet2Vol;
  G4VPhysicalVolume* FingerCounter;
  G4VPhysicalVolume* AnaBar;

  DetectorSD*        fDetSD;
  PMTSD*             fPMTSD;

  G4int              fTumourOn;
  G4double           fTumourRadius;
  G4double           fTumourHeight;
  G4double           fAnaBarXpos;

};
#endif

//---------------------------------------------------------------------------

