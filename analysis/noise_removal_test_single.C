#include <iostream>
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "AddNoise.C"
#include "ChiSquareTest.C"

using namespace std;

void noise_removal_test_single() {

    gSystem->Exec("cp data/AnaBarMC_3007_orig.root data/AnaBarMC_3007.root");
    AddNoise(3007, 1);
    ChiSquareTest(3007);

    return;
}

