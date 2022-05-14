
class Features {
public:
    Int_t iEvt;            // event index
    Float_t q2;            // lepton momentum sqaured
    Float_t miss2;         // missing mass squared
    Float_t pB;            // b-hadron momentum
    Float_t EB;            // b-hadron energy
    Float_t pHc;           // c-hadron momentum
    Float_t EHc;           // c-hadron energy
    Float_t pMu;           // muon momentum
    Float_t EMu;           // muon energy
    Float_t sMinMuBVert;   // closest distance between Hc track and muon track (=between deduced B edcay vertex and muon track)
    Float_t sMinMuHcVert;  // closest distance between Hc decay vertex and muon track
    Float_t sMinMuTr;      // closest track to muon
    Float_t sMinHcTr;      // closest track to reconstructed Hc
    Float_t sPVHc;         // distance from PV to Hc decay vertex
    Float_t mHcMu;         // invariant mass of Hc + mu
    Float_t pPerp;
    // Float_t pPerpHc;
    Float_t mCorr;
    Float_t D0Max;
    Float_t DzMax;
    Float_t D0Sum;
    Float_t DzSum;
    Float_t ENeutral03;
    Float_t ENeutral06;
    Float_t ENeutral03Hadron;
    Float_t ENeutral06Hadron;
    Float_t ENeutral03Photon;
    Float_t ENeutral06Photon;
    Float_t ECharge03;
    Float_t ECharge06;
    Float_t ECharge03PV;
    Float_t ECharge06PV;
    Float_t ECharge03DV;
    Float_t ECharge06DV;
    Float_t mK0SHcMu;  // invariant mass of the reconstructed K0S + Hc + mu;
    Float_t pK0S;      // momentum of the reconstructed K0S

    Float_t q2True;
    Float_t miss2True;
    Float_t EBTrue;
    Float_t pBTrue;
    Float_t sMinMuHcVertTrue;
};
