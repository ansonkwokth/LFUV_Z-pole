
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
    Float_t pPerp;         // from arXiv:2001.03225v2 eqn3
    // Float_t pPerpHc;
    Float_t mCorr;             // from arXiv:2001.03225v2 eqn3
    Float_t D0Max;             // max. of particle transverse impact parameter
    Float_t DzMax;             // max. of particle longitudinal impact parameter
    Float_t D0Sum;             // sum of particle transverse impact parameter
    Float_t DzSum;             // sum of particle longitudinal impact parameter
    Float_t ENeutral03;        // sum of neutral particle energ, which has angle between B < 0.3rad
    Float_t ENeutral06;        // sum of neutral particle energ, which has angle  between B < 0.6rad
    Float_t ENeutral03Hadron;  // sum of neutral hadron energ, which has angle between B < 0.3rad
    Float_t ENeutral06Hadron;  // sum of neutral hadron energ, which has angle between B < 0.6rad
    Float_t ENeutral03Photon;  // sum of photon energ, which has angle between B < 0.6rad
    Float_t ENeutral06Photon;  // sum of photon energ, which has angle between B < 0.6rad
    Float_t ECharge03;         // sum of charged  particle energ, which has angle between B < 0.6rad
    Float_t ECharge06;         // sum of charged particle energ, which has angle between B < 0.6rad
    Float_t ECharge03PV;       // sum of chaged PV particle energ, which has angle between B < 0.6rad
    Float_t ECharge06PV;       // sum of charged PV particle energ, which has angle between B < 0.6rad
    Float_t ECharge03DV;       // sum of charged DV particle energ, which has angle between B < 0.6rad
    Float_t ECharge06DV;       // sum of charged DV particle energ, which has angle between B < 0.6rad
    Float_t mK0SHcMu;          // invariant mass of the reconstructed K0S + Hc + mu;
    Float_t pK0S;              // momentum of the reconstructed K0S
    Float_t DeltaM;            // m(Ds gamma) - m(Ds)
    Float_t EPho;
    Float_t etaPho;
    Float_t phiPho;

    Float_t q2True;
    Float_t miss2True;
    Float_t EBTrue;
    Float_t pBTrue;
    Float_t sMinMuHcVertTrue;
    Float_t sMinMuBVertTrue;
    Int_t correctPhoton;
    Float_t EPhoTrue;  // truth photon energy
    Float_t etaPhoTrue;
    Float_t phiPhoTrue;
    Int_t isDsPho;            // check if the cheated photon is the Ds photon
    Float_t DeltaRDsPhoTrue;  // angular separation between truth Ds and Photon
    Int_t bkgFromDsstar;

    Float_t DeltaRPhoPhoTrue;  // angular separation of the true and reco photon
                               //

    Float_t ETower;
};
