#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"

//---------------------------------------------------------------------------

class WLSMaterials;
class G4Material;

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
  inline G4VPhysicalVolume* GetDet3Vol()    { return fDet3Vol;  };
  inline G4VPhysicalVolume* GetDet15Vol()    { return fDet15Vol;  };
  inline G4VPhysicalVolume* GetFingerVol()    { return FingerCounter;  };
  inline G4VPhysicalVolume* GetAnaBar1Vol()    { return AnaBar1;  };
  inline G4VPhysicalVolume* GetAnaBar2Vol()    { return AnaBar2;  };
  inline G4VPhysicalVolume* GetAnaBar3Vol()    { return AnaBar3;  };
  inline G4VPhysicalVolume* GetClad1Vol()    { return physiClad1;  };
  inline G4VPhysicalVolume* GetClad2Vol()    { return physiClad2;  };
  inline G4VPhysicalVolume* GetClad3Vol()    { return physiClad3;  };
  inline G4VPhysicalVolume* GetWLSfiber1Vol()    { return physiWLSfiber1;  };
  inline G4VPhysicalVolume* GetWLSfiber2Vol()    { return physiWLSfiber2;  };
  inline G4VPhysicalVolume* GetWLSfiber3Vol()    { return physiWLSfiber3;  };
  inline DetectorSD*        GetDetSD()      { return fDetSD;    };
  inline PMTSD*             GetPMTSD()      { return fPMTSD;    };

  inline void SetTumourOn     ( G4int    tumon )   { fTumourOn     = tumon; }
  inline void SetTumourRadius ( G4double radius )  { fTumourRadius = radius; }
  inline void SetTumourHeight ( G4double height )  { fTumourHeight = height; }
  inline void SetAnaBarXpos   ( G4double AnaBarXpos )    { fAnaBarXpos   = AnaBarXpos; }
  inline void SetAnaBarLength   ( G4double AnaBarLength )    { fAnaBarLength   = AnaBarLength; }
  inline void SetAnaBarWidth   ( G4double AnaBarWidth )    { fAnaBarWidth   = AnaBarWidth; }
  inline void SetAnaBarThickness   ( G4double AnaBarThickness )    { fAnaBarThickness   = AnaBarThickness; }
  inline void SetFingerLength   ( G4double FingerLength )    { fFingerLength   = FingerLength; }
  inline void SetFingerWidth   ( G4double FingerWidth )    { fFingerWidth   = FingerWidth; }
  inline void SetFingerThickness   ( G4double FingerThickness )    { fFingerThickness   = FingerThickness; }
  inline void SetFibreDiameter   ( G4double FibreDiameter )    { fFibreDiameter  = FibreDiameter; }
  inline void SetFibreLength   ( G4double FibreLength )    { fFibreLength  = FibreLength; }
  inline void SetHoleDiameter   ( G4double HoleDiameter )    { fHoleDiameter  = HoleDiameter; }
  inline void SetHoleLength   ( G4double HoleLength )    { fHoleLength  = HoleLength; }
  inline void SetCladdingDiameter   ( G4double CladdingDiameter )    { fCladdingDiameter  = CladdingDiameter; }
  inline void SetCladdingLength   ( G4double CladdingLength )    { fCladdingLength  = CladdingLength; }
  inline void SetPhotoCathodeDiameter   ( G4double PhotoCathodeDiameter )    { fPhotoCathodeDiameter  = PhotoCathodeDiameter; }
  inline void SetPhotoCathodeThickness   ( G4double PhotoCathodeThickness)    { fPhotoCathodeThickness  = PhotoCathodeThickness; }

  G4Material* FindMaterial(G4String);
  
  private:

  WLSMaterials* fMaterials;

  G4NistManager*     fNistManager;
  DetectorMessenger* fDetMessenger;

  G4VPhysicalVolume* fExpHall;
  G4VPhysicalVolume* fDet1Vol;
  G4VPhysicalVolume* fDet2Vol;
  G4VPhysicalVolume* fDet3Vol;
  G4VPhysicalVolume* fDet15Vol;
  G4VPhysicalVolume* FingerCounter;
  G4VPhysicalVolume* AnaBar1;
  G4VPhysicalVolume* AnaBar2;
  G4VPhysicalVolume* AnaBar3;
  G4VPhysicalVolume* physiClad1;
  G4VPhysicalVolume* physiClad2;
  G4VPhysicalVolume* physiClad3;
  G4VPhysicalVolume* physiWLSfiber1;
  G4VPhysicalVolume* physiWLSfiber2;
  G4VPhysicalVolume* physiWLSfiber3;

  DetectorSD*        fDetSD;
  PMTSD*             fPMTSD;

  G4int              fTumourOn;
  G4double           fTumourRadius;
  G4double           fTumourHeight;
  G4double           fAnaBarXpos;

  G4double           fAnaBarLength;
  G4double           fAnaBarWidth;
  G4double           fAnaBarThickness;
  G4double           fFingerLength;
  G4double           fFingerWidth;
  G4double           fFingerThickness;
  G4double           fFibreDiameter;
  G4double           fFibreLength;
  G4double           fHoleDiameter;
  G4double           fHoleLength;
  G4double           fCladdingDiameter;
  G4double           fCladdingLength;
  G4double           fPhotoCathodeDiameter;
  G4double           fPhotoCathodeThickness;

};
#endif

//---------------------------------------------------------------------------

