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
  inline G4VPhysicalVolume* GetDet1Vol_rm()    { return fDet1Vol_rm;  };
  inline G4VPhysicalVolume* GetDet1Vol_lm()    { return fDet1Vol_lm;  };
  inline G4VPhysicalVolume* GetDet1Vol_mu()    { return fDet1Vol_mu;  };
  inline G4VPhysicalVolume* GetDet1Vol_ml()    { return fDet1Vol_ml;  };
  inline G4VPhysicalVolume* GetDet1Vol_ru()    { return fDet1Vol_ru;  };
  inline G4VPhysicalVolume* GetDet1Vol_lu()    { return fDet1Vol_lu;  };
  inline G4VPhysicalVolume* GetDet1Vol_rl()    { return fDet1Vol_rl;  };
  inline G4VPhysicalVolume* GetDet1Vol_ll()    { return fDet1Vol_ll;  };
  inline G4VPhysicalVolume* GetDet2Vol()    { return fDet2Vol;  };
  inline G4VPhysicalVolume* GetFingerVol()    { return FingerCounter;  };
  inline G4VPhysicalVolume* GetAnaBarVol()    { return AnaBar;  };
  inline G4VPhysicalVolume* GetAnaBarVol_rm()    { return AnaBar_rm;  };
  inline G4VPhysicalVolume* GetAnaBarVol_lm()    { return AnaBar_lm;  };
  inline G4VPhysicalVolume* GetAnaBarVol_mu()    { return AnaBar_mu;  };
  inline G4VPhysicalVolume* GetAnaBarVol_ml()    { return AnaBar_ml;  };
  inline G4VPhysicalVolume* GetAnaBarVol_ru()    { return AnaBar_ru;  };
  inline G4VPhysicalVolume* GetAnaBarVol_lu()    { return AnaBar_lu;  };
  inline G4VPhysicalVolume* GetAnaBarVol_rl()    { return AnaBar_rl;  };
  inline G4VPhysicalVolume* GetAnaBarVol_ll()    { return AnaBar_ll;  };
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
  G4VPhysicalVolume* fDet1Vol_rm;
  G4VPhysicalVolume* fDet1Vol_lm;
  G4VPhysicalVolume* fDet1Vol_mu;
  G4VPhysicalVolume* fDet1Vol_ml;
  G4VPhysicalVolume* fDet1Vol_ru;
  G4VPhysicalVolume* fDet1Vol_lu;
  G4VPhysicalVolume* fDet1Vol_rl;
  G4VPhysicalVolume* fDet1Vol_ll;
  G4VPhysicalVolume* fDet2Vol;
  G4VPhysicalVolume* FingerCounter;
  G4VPhysicalVolume* AnaBar;
  G4VPhysicalVolume* AnaBar_rm;
  G4VPhysicalVolume* AnaBar_lm;
  G4VPhysicalVolume* AnaBar_mu;
  G4VPhysicalVolume* AnaBar_ml;
  G4VPhysicalVolume* AnaBar_ru;
  G4VPhysicalVolume* AnaBar_lu;
  G4VPhysicalVolume* AnaBar_rl;
  G4VPhysicalVolume* AnaBar_ll;

  DetectorSD*        fDetSD;
  PMTSD*             fPMTSD;

  G4int              fTumourOn;
  G4double           fTumourRadius;
  G4double           fTumourHeight;
  G4double           fAnaBarXpos;

};
#endif

//---------------------------------------------------------------------------

