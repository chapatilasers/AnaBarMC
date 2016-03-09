#ifndef SteppingAction_H
#define SteppingACtion_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

//---------------------------------------------------------------------------

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  ~SteppingAction();
  virtual void UserSteppingAction(const G4Step*);

  inline void SetOneStepPrimaries(G4bool b) { fOneStepPrimaries = b;    }
  inline G4bool GetOneStepPrimaries()       { return fOneStepPrimaries; }
  
private:
  G4bool fOneStepPrimaries;

};

#endif

//---------------------------------------------------------------------------
