#pragma once
#include "Rtypes.h"

struct iFinalStates {
    Int_t foundFromC = 0;  // check if all decay products from c-hadron are found
    Int_t foundAll = 0;    // check if all decay productes from b-hadron are found
    Int_t iP = 99999;      // index of proton
    Int_t iK = 99999;      // index of kaon
    Int_t iPi = 99999;     // index of pion
    Int_t iMu = 99999;     // index of muon
};

