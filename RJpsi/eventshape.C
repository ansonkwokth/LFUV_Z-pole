#include <cmath>
#include <iostream>
using namespace std;

#include "Rtypes.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

R__ADD_LIBRARY_PATH($DELPHES)
R__LOAD_LIBRARY(libDelphes)

// find the final states of the decay{{{
struct iFinalStates {
    Int_t foundFromC = 0;
    Int_t foundAll = 0;
    Int_t iMuPos = 99999;
    Int_t iMuNeg = 99999;
    Int_t iMu = 99999;
};
//}}}

// return if the input is signal {{{
Int_t ClassifySignal(TClonesArray* branchParticle, TLorentzVector* BTrue, TLorentzVector* HcTrue, TLorentzVector* muTrue, TVector3* v3HcTrue, TVector3* v3MuTrue, iFinalStates* iFSTrue, Int_t mode) {
    GenParticle* particle;
    GenParticle* particleB;
    GenParticle* particleC;
    GenParticle* particle0;
    GenParticle* particle0M;
    GenParticle* particle0MM;
    Int_t nParticles = branchParticle->GetEntries();

    Int_t iBTrue, foundBTrue = 0;
    Int_t iCTrue, foundCTrue = 0;
    Int_t iC1True, foundC1True = 0;
    Int_t iC2True, foundC2True = 0;
    Int_t iMuTrue, foundMuTrue = 0;
    for (int ipB = 0; ipB < nParticles; ipB++) {
        foundBTrue = 0;
        particleB = (GenParticle*)branchParticle->At(ipB);
        if (abs(particleB->PID) == 541) {
            iBTrue = ipB;
            foundBTrue = 1;

            Int_t iCTrue_i, foundCTrue_i = 0;
            Int_t iC1True_i, foundC1True_i = 0;
            Int_t iC2True_i, foundC2True_i = 0;
            Int_t iMuTrue_i, foundMuTrue_i = 0;
            for (int ip2 = 0; ip2 < nParticles; ip2++) {
                particle0 = (GenParticle*)branchParticle->At(ip2);
                if (particle0->M1 == -1) {
                    continue;
                }
                particle0M = (GenParticle*)branchParticle->At(particle0->M1);
                if (particle0M->M1 == -1) {
                    continue;
                }
                particle0MM = (GenParticle*)branchParticle->At(particle0M->M1);

                if (abs(particle0->PID) == 443 && particle0->M1 == iBTrue) {
                    iCTrue_i = ip2;
                    foundCTrue_i = 1;
                }
                if (particle0->PID == -13 && abs(particle0M->PID) == 443 && particle0M->M1 == iBTrue) {
                    iC1True_i = ip2;
                    foundC1True_i = 1;
                }
                if (particle0->PID == 13 && abs(particle0M->PID) == 443 && particle0M->M1 == iBTrue) {
                    iC2True_i = ip2;
                    foundC2True_i = 1;
                }

                if (mode == 1 && (abs(particle0->PID) == 13 && abs(particle0M->PID) == 15 && particle0M->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
                if (mode == 2 && (abs(particle0->PID) == 13 && particle0->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
            }
            // cout <<  foundCTrue_i <<"; "<< foundC1True_i <<"; "<<
            // foundC2True_i << "; " << foundMuTrue_i<<endl;
            if (foundCTrue_i == 1 && foundC1True_i == 1 && foundC2True_i == 1 && foundMuTrue_i == 1) {
                iCTrue = iCTrue_i;
                foundCTrue = 1;
                iC1True = iC1True_i;
                foundC1True = 1;
                iC2True = iC2True_i;
                foundC2True = 1;
                iMuTrue = iMuTrue_i;
                foundMuTrue = 1;
                break;
            }
        }
    }
    // cout <<  foundBTrue <<"; "<< foundCTrue <<"; "<< foundCTrue <<"; "<<
    // foundC2True <<"; "<< foundC3True << "; " << foundMuTrue <<endl;
    if (foundBTrue == 1 && foundCTrue == 1 && foundCTrue == 1 && foundC2True == 1 && foundMuTrue == 1) {
        iFinalStates iFSTrue_;
        iFSTrue_.iMuPos = iC1True;
        iFSTrue_.iMuNeg = iC2True;
        iFSTrue_.iMu = iMuTrue;

        TLorentzVector BTrue_, HcTrue_, muTrue_;
        GenParticle* particleHc;
        GenParticle* particleMu;
        GenParticle* particleC1;
        particleB = (GenParticle*)branchParticle->At(iBTrue);
        BTrue_.SetPtEtaPhiE(particleB->PT, particleB->Eta, particleB->Phi, particleB->E);
        particleHc = (GenParticle*)branchParticle->At(iCTrue);
        HcTrue_.SetPtEtaPhiE(particleHc->PT, particleHc->Eta, particleHc->Phi, particleHc->E);
        particleMu = (GenParticle*)branchParticle->At(iMuTrue);
        muTrue_.SetPtEtaPhiE(particleMu->PT, particleMu->Eta, particleMu->Phi, particleMu->E);

        particleC1 = (GenParticle*)branchParticle->At(iC1True);
        TVector3 v3HcTrue_, v3MuTrue_;
        v3HcTrue_.SetXYZ(particleC1->X, particleC1->Y, particleC1->Z);
        v3MuTrue_.SetXYZ(particleMu->X, particleMu->Y, particleMu->Z);

        *iFSTrue = iFSTrue_;
        *BTrue = BTrue_;
        *HcTrue = HcTrue_;
        *muTrue = muTrue_;
        *v3HcTrue = v3HcTrue_;
        *v3MuTrue = v3MuTrue_;

        return 1;
    } else {
        return 0;
    }
}
// }}}

// Return if the event is the targeted bkg (from inclusive samples) {{{
Int_t ClassifyBkg(TClonesArray* branchParticle, const string type) {
    // Inclusive samepls (for Comb+Cascade, Comb, Casde, Inclusive)
    GenParticle* particle;
    Int_t nParticles = branchParticle->GetEntries();

    //=== check c-hadron ===
    Int_t nPosHc = 0;  // number of positive targeted c-hadron in the truth
                       // level (Jpsi, Ds, Lambdac)
    for (int ip = 0; ip < nParticles; ip++) {
        particle = (GenParticle*)branchParticle->At(ip);
        if (abs(particle->PID) == 443) {
            nPosHc += 1;
        }  // check number of pos lambdac
    }
    if (not(nPosHc >= 1)) {
        return 0;
    }  // must have at least one lambdac

    //=== check muons ===
    Int_t nPosMu = 0;  // number of positive muons
    Int_t nNegMu = 0;  // number of negative muons
    for (int ip = 0; ip < nParticles; ip++) {
        particle = (GenParticle*)branchParticle->At(ip);
        if (particle->PID == -13) {
            nPosMu += 1;
        }
        if (particle->PID == 13) {
            nNegMu += 1;
        }
    }
    // total must be at least 3, (2 from Jpsi, 1 unapired), and at least 1 pos,
    // 1 neg from Jpsi
    if (not(nPosMu > 0 && nNegMu > 0 && nPosMu + nNegMu >= 3)) {
        return 0;
    }

    //=== check the decay chain ===
    Int_t isLeptonFromB = 0;      // check if the lepton is from b-hadron
    Int_t isLeptonSemiFromB = 0;  // check if the lepton is semileptonically decay from b-hadron
    Int_t nLeptonStepFromB = 0;
    Int_t isCFromB = 0;  // check if c-hadron is from B
    Int_t isSignal = 0;  // chekc if event has signal in it

    GenParticle* particle1;
    GenParticle* particle1M;  // particle1's mother
    GenParticle* particle2;
    GenParticle* particle2M;  // particle2's mother
    GenParticle* particleL;   // lepton
    GenParticle* particleC;   // c-hadron
    for (int ip1 = 0; ip1 < nParticles; ip1++) {
        particle1 = (GenParticle*)branchParticle->At(ip1);
        // identify the muon
        if (abs(particle1->PID) == 13) {
            isLeptonFromB = 0;
            isLeptonSemiFromB = 0;
            isCFromB = 0;

            // mother of particle0
            particle1M = (GenParticle*)branchParticle->At(particle1->M1);
            if (abs(particle1M->PID) == 443) {
                continue;
            }

            // identify the lepton (mu/tau)
            Int_t iLepton = ip1;
            if (abs(particle1M->PID) == 15) {
                iLepton = particle1->M1;
            }
            particleL = (GenParticle*)branchParticle->At(iLepton);
            Int_t iLeptonMother = particleL->M1;
            particle1M = (GenParticle*)branchParticle->At(iLeptonMother);
            // cout << "Lepton: " << particleL->PID << " (" << iLepton << ")
            // from: ";

            // finding the b-hadron that lepton from
            Int_t BHadron_idx = 99999;
            nLeptonStepFromB = 0;
            while (true) {
                // cout << particle1M->PID << " < ";
                if (int(abs(particle1M->PID) / 100) == 5 || int(abs(particle1M->PID) / 1000) == 5 || int(abs(particle1M->PID) / 10000) == 5 || int((abs(particle1M->PID) % 1000) / 100) == 5) {
                    BHadron_idx = iLeptonMother;
                    // cout << "\n From b-hadron (IDX): " << particle1M->PID <<
                    // " (" << BHadron_idx << ")" << endl;
                    break;
                } else if (int(particle1M->PID / 10) == 0) {
                    // cout << "\nNo b-hadron, terminated at: " <<
                    // particle1M->PID << endl;
                    break;
                }
                nLeptonStepFromB += 1;
                iLeptonMother = particle1M->M1;
                particle1M = (GenParticle*)branchParticle->At(iLeptonMother);
            }
            if (BHadron_idx != 99999) {
                isLeptonFromB = 1;
            }
            if (BHadron_idx != 99999 && nLeptonStepFromB == 0) {
                isLeptonSemiFromB = 1;
            }

            // if lepton is from b-hadron, then find if c-hadron is from the
            // same mother
            if (isLeptonFromB == 1) {
                for (int ip2 = 0; ip2 < nParticles; ip2++) {
                    particle2 = (GenParticle*)branchParticle->At(ip2);
                    // cout<< ip2<<" ;"<<abs(particle2->PID) <<endl;
                    if (abs(particle2->PID) == 443) {
                        Int_t iC = ip2;
                        particleC = (GenParticle*)branchParticle->At(iC);
                        Int_t iCMother = particleC->M1;
                        particle2M = (GenParticle*)branchParticle->At(iCMother);

                        if (iCMother == BHadron_idx && isLeptonSemiFromB == 1) {
                            isSignal = 1;
                            continue;
                        }
                        // cout << "c-hadron " << particleC->PID << " (" <<
                        // iCMother << ") from: ";

                        while (true) {
                            // cout << particle2M->PID << " (" << iCMother <<
                            // ")" << " < ";
                            if (iCMother == BHadron_idx) {
                                // cout << "/nc-hadron from the smae b-hadron: "
                                // << particle2M->PID << endl;
                                isCFromB = 1;
                                break;
                            }
                            if (int(particle2M->PID / 10) == 0) {
                                // cout << "\nc-hadron not from same b-hadron
                                // terminate at: " << particle2M->PID << endl;
                                break;
                            }
                            iCMother = particle2M->M1;
                            particle2M = (GenParticle*)branchParticle->At(iCMother);
                        }
                        if (isCFromB == 1 || isSignal == 1) {
                            break;
                        }
                    }
                }
                if (isCFromB == 1 || isSignal == 1) {
                    break;
                }
            }
        }
    }
    if (isSignal == 1) {
        return 0;
    }

    Int_t isComb = 0;
    Int_t isCascade = 0;
    Int_t isInclusive = 0;

    if (isLeptonSemiFromB == 1 && isCFromB == 1) {
        isInclusive = 1;
    }
    if (isLeptonFromB == 1 && isLeptonSemiFromB == 0 && isCFromB == 1) {
        isCascade = 1;
    }
    if (isCFromB == 0) {
        isComb = 1;
    }

    if (isInclusive + isCascade + isComb != 1) {
        cout << "HAVE BUG IN CLASSIFYING BKG!";
    }

    // cout<<  "" <<endl;
    // cout << "Signal: " << isSignal << "; Comb:" << isComb  << "; Cascade: "
    // << isCascade<< "; Inclsuive: " << isInclusive << endl;
    if (type == "b1") {
        if (isCascade == 0 && isComb == 0) { /*cout<<  "???" <<endl;*/
            return 0;
        } else {
            return 1;
        }
    }
    if (type == "b2") {
        if (isComb == 0) {
            return 0;
        } else {
            return 1;
        }
    }
    if (type == "b3") {
        if (isCascade == 0) {
            return 0;
        } else {
            return 1;
        }
    }
    if (type == "b4") {
        if (isInclusive == 0) {
            return 0;
        } else {
            return 1;
        }
    } else {
        return 0;
    }
}

///}}}

// Return if the event is misID bkg // {{{
Int_t ClassifyMisID(TClonesArray* branchParticle) {
    // skip signal events
    TLorentzVector dummy;
    TVector3 dummy2;
    iFinalStates dummy3;
    if (ClassifySignal(branchParticle, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy3, 1) == 1) {
        return 0;
    } else if (ClassifySignal(branchParticle, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy3, 2) == 1) {
        return 0;
    } else {
        GenParticle* particle;
        Int_t nParticles = branchParticle->GetEntries();

        //=== check c-hadron ===
        Int_t nHc = 0;  // number of targeted c-hadron in the truth level (Jpsi,
                        // Ds, Lambdac)
        for (int ip = 0; ip < nParticles; ip++) {
            particle = (GenParticle*)branchParticle->At(ip);
            if (abs(particle->PID) == 443) {
                nHc += 1;
            }  // check number of Jpsi
        }
        if (not(nHc >= 1)) {
            return 0;
        }  // must have at least one Jpsi

        //=== check pions ===
        Int_t nPi = 0;  // number of pions (displaced)
        for (int ip = 0; ip < nParticles; ip++) {
            particle = (GenParticle*)branchParticle->At(ip);
            if (abs(particle->PID) == 211 && (not(particle->X == 0 && particle->Y == 0 && particle->Z == 0))) {
                nPi += 1;
            }
        }
        if (not(nPi >= 1)) {
            return 0;
        }
        return 1;
    }
}

//}}}

//{{{
Float_t calCosinceTheta(TLorentzVector p1, TLorentzVector p2) {
    return (p1.Px() * p2.Px() + p1.Py() * p2.Py() + p1.Pz() * p2.Pz()) / (p1.P() * p2.P());
}
//}}}

class FWMoments {
public:
    Int_t iEvt;
    Float_t H_EE0 = 0;
    Float_t H_EE1 = 0;
    Float_t H_EE2 = 0;
    Float_t H_EE3 = 0;
    Float_t H_EE4 = 0;
    Float_t H_EE5 = 0;
    Float_t H_EE6 = 0;
    Float_t H_EE7 = 0;
    Float_t H_EE8 = 0;
    Float_t H_EE9 = 0;
    Float_t H_EE10 = 0;
};

FWMoments updateFWMs(Float_t E1, Float_t E2, Float_t cosine_theta, FWMoments FWMs) {
    FWMoments FWMs_;
    FWMs_.H_EE0 = FWMs.H_EE0 + E1 * E2;
    FWMs_.H_EE1 = FWMs.H_EE1 + E1 * E2 * cosine_theta;
    FWMs_.H_EE2 = FWMs.H_EE2 + E1 * E2 * (0.5 * (3 * pow(cosine_theta, 2) - 1));
    FWMs_.H_EE3 = FWMs.H_EE3 + E1 * E2 * (0.5 * (5 * pow(cosine_theta, 3) - 3 * cosine_theta));
    FWMs_.H_EE4 = FWMs.H_EE4 + E1 * E2 * ((1.0 / 8.0) * (35 * pow(cosine_theta, 4) - 30 * pow(cosine_theta, 2) + 3));
    FWMs_.H_EE5 = FWMs.H_EE5 + E1 * E2 * ((1.0 / 8.0) * (63 * pow(cosine_theta, 5) - 70 * pow(cosine_theta, 3) + 15 * cosine_theta));
    FWMs_.H_EE6 = FWMs.H_EE6 + E1 * E2 * ((1.0 / 16.0) * (231 * pow(cosine_theta, 6) - 315 * pow(cosine_theta, 4) + 105 * pow(cosine_theta, 2) - 5));
    FWMs_.H_EE7 = FWMs.H_EE7 + E1 * E2 * ((1.0 / 16.0) * (429 * pow(cosine_theta, 7) - 693 * pow(cosine_theta, 5) + 315 * pow(cosine_theta, 3) - 35 * cosine_theta));
    FWMs_.H_EE8 = FWMs.H_EE8 + E1 * E2 * ((1.0 / 128.0) * (6435 * pow(cosine_theta, 8) - 12012 * pow(cosine_theta, 6) + 6930 * pow(cosine_theta, 4) - 1260 * pow(cosine_theta, 2) + 35));
    FWMs_.H_EE9 = FWMs.H_EE9 + E1 * E2 * ((1.0 / 128.0) * (12155 * pow(cosine_theta, 9) - 25740 * pow(cosine_theta, 7) + 18018 * pow(cosine_theta, 5) - 4620 * pow(cosine_theta, 3) + 315 * cosine_theta));
    FWMs_.H_EE10 = FWMs.H_EE10 + E1 * E2 * ((1.0 / 256) * (46189 * pow(cosine_theta, 10) - 109395 * pow(cosine_theta, 8) + 90090 * pow(cosine_theta, 6) - 30030 * pow(cosine_theta, 4) + 3465 * pow(cosine_theta, 2) - 63));

    return FWMs_;
}

Int_t nEvt = 0;
void eventshape(const string type, const Bool_t save = true, Int_t num_test = 0) {
    cout << "\n\n\n\n\n\n\n\n\n\n";

    string typeName;
    const char* inputFile;
    const char* outputFile;
    if (type == "s1") {
        typeName = "Bc->Jpsi tau nu. ";
        cout << "Bc->Jpsi tau nu. " << endl;
        inputFile = "~/Projects/LFUV/RJpsi/Bcjpsitaunu_50m.root";
        outputFile = "~/Projects/LFUV/RJpsi/features/JpsiTauNu_FWM.root";
    } else if (type == "s2") {
        cout << "Bc->Jpsi mu nu. " << endl;
        inputFile = "~/Projects/LFUV/RJpsi/BcJpsimunu0-2.root";
        outputFile = "~/Projects/LFUV/RJpsi/features/JpsiMuNu_FWM.root";
    } else if (type == "b1") {
        cout << "Comb+Cascade Bkg. " << endl;
        inputFile = "~/Projects/LFUV/RJpsi/RJpsi_comb_200m_seed1.root";
        outputFile = "~/Projects/LFUV/RJpsi/features/RJpsiCombCascade_FWM.root";
    } else if (type == "b2") {
        cout << "Comb Bkg. " << endl;
        inputFile = "~/Projects/LFUV/RJpsi/RJpsi_comb_200m_seed1.root";
        outputFile = "~/Projects/LFUV/RJpsi/features/RJpsiComb_FWM.root";
    } else if (type == "b3") {
        cout << "Cascade Bkg. " << endl;
        inputFile = "~/Projects/LFUV/RJpsi/RJpsi_comb_200m_seed1.root";
        outputFile = "~/Projects/LFUV/RJpsi/features/RJpsiCascade_FWM.root";
    } else if (type == "b4") {
        cout << "Inclusive Bkg. " << endl;
        inputFile = "~/Projects/LFUV/RJpsi/RJpsi_comb_200m_seed1.root";
        outputFile = "~/Projects/LFUV/RJpsi/features/RJpsiInclusive_FWM.root";
    } else if (type == "b5") {
        cout << "MisID Bkg. " << endl;
        inputFile = "~/Projects/LFUV/RJpsi/RJpsi_comb_200m_seed1.root";
        outputFile = "~/Projects/LFUV/RJpsi/features/RJpsiMisID_FWM.root";
    } else {
        cout << "Not match. ";
    }

    // Load lib, and read data
    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    chain.Add(inputFile);
    ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
    Int_t numberOfEntries = treeReader->GetEntries();

    // clone the branches
    TClonesArray* branchParticle = treeReader->UseBranch("Particle");
    TClonesArray* branchTrack = treeReader->UseBranch("Track");
    TClonesArray* branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray* branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

    // defining constants
    const Float_t mK0 = 0.497661;

    GenParticle* particle;
    Track* track;
    Tower* eflowphoton;
    Tower* eflowneutralhadron;

    if (not save) outputFile = "dummy.root";

    TFile fea(outputFile, "recreate");
    TTree tr("t", "features");
    FWMoments* FWMs = new FWMoments;
    // Features features;
    tr.Branch("FWMs", &FWMs);

    cout << "\nReconstructing: " << typeName << endl;
    cout << "Reading from:\t" << inputFile << endl;
    cout << "Save: " << save << endl;
    cout << "Writing to:\t" << outputFile << endl;
    cout << "# of events:\t" << numberOfEntries << "\n\n";

    if (num_test == 0) num_test = numberOfEntries;
    for (Int_t i_en = 0; i_en < numberOfEntries; i_en++) {
        if ((i_en % 10000) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << "\n";
        if (i_en >= num_test) break;

        treeReader->ReadEntry(i_en);  // reading the entry

        //==========================================================================
        //===============   Classifying in event type in truth level ===============
        //==========================================================================
        iFinalStates iFSTrue;
        TLorentzVector BTrue, HcTrue, muTrue;
        TVector3 v3HcTrue, v3MuTrue;
        Int_t passing = 0;

        if (type == "s1") {
            passing = ClassifySignal(branchParticle, &BTrue, &HcTrue, &muTrue, &v3HcTrue, &v3MuTrue, &iFSTrue, 1);
        } else if (type == "s2") {
            passing = ClassifySignal(branchParticle, &BTrue, &HcTrue, &muTrue, &v3HcTrue, &v3MuTrue, &iFSTrue, 2);
        } else if (type == "b1" || type == "b1" || type == "b2" || type == "b3" || type == "b4") {
            passing = ClassifyBkg(branchParticle, type);
        } else if (type == "b5") {
            passing = ClassifyMisID(branchParticle);
        } else {
            cout << " No such Signal/Bkg" << endl;
        }
        if (passing == 0) continue;

        Int_t nTracks = branchTrack->GetEntries();
        Track* track1;
        Track* track2;
        Int_t nPhotons = branchEFlowPhoton->GetEntries();
        Tower* photon1;
        Tower* photon2;
        Int_t nNeutral = branchEFlowNeutralHadron->GetEntries();
        Tower* neutral1;
        Tower* neutral2;

        nEvt += 1;
        // cout << "---------------------Event: " << i_en << endl;

        FWMoments FWMs_;
        // 1st tracks loop
        for (Int_t itr1 = 0; itr1 < nTracks; itr1++) {
            track1 = (Track*)branchTrack->At(itr1);
            TLorentzVector cp1;
            cp1.SetPtEtaPhiE(track1->PT, track1->Eta, track1->Phi, track1->P);
            Float_t E1 = cp1.E();

            // tracks loop
            for (Int_t itr2 = 0; itr2 < nTracks; itr2++) {
                track2 = (Track*)branchTrack->At(itr2);
                TLorentzVector cp2;
                cp2.SetPtEtaPhiE(track2->PT, track2->Eta, track2->Phi, track2->P);
                Float_t E2 = cp2.E();
                Float_t cosine_theta = calCosinceTheta(cp1, cp2);
                FWMs_ = updateFWMs(E1, E2, cosine_theta, FWMs_);
            }
            // neutral hadron loop
            for (Int_t inh2 = 0; inh2 < nNeutral; inh2++) {
                neutral2 = (Tower*)branchEFlowNeutralHadron->At(inh2);
                Float_t pt2 = neutral2->ET * pow(neutral2->E * neutral2->E - mK0 * mK0, 0.5) / neutral2->E;
                TLorentzVector nh2;
                nh2.SetPtEtaPhiE(pt2, neutral2->Eta, neutral2->Phi, neutral2->E);
                Float_t E2 = nh2.E();
                Float_t cosine_theta = calCosinceTheta(cp1, nh2);
                FWMs_ = updateFWMs(E1, E2, cosine_theta, FWMs_);
            }
            // photon loop
            for (Int_t iph2 = 0; iph2 < nPhotons; iph2++) {
                photon2 = (Tower*)branchEFlowPhoton->At(iph2);
                TLorentzVector ph2;
                ph2.SetPtEtaPhiE(photon2->ET, photon2->Eta, photon2->Phi, photon2->E);
                Float_t E2 = ph2.E();
                Float_t cosine_theta = calCosinceTheta(cp1, ph2);
                FWMs_ = updateFWMs(E1, E2, cosine_theta, FWMs_);
            }
        }

        // 1st neutral hadron loop
        for (Int_t inh1 = 0; inh1 < nNeutral; inh1++) {
            neutral1 = (Tower*)branchEFlowNeutralHadron->At(inh1);
            Float_t pt1 = neutral1->ET * pow(neutral1->E * neutral1->E - mK0 * mK0, 0.5) / neutral1->E;
            TLorentzVector nh1;
            nh1.SetPtEtaPhiE(pt1, neutral1->Eta, neutral1->Phi, neutral1->E);
            Float_t E1 = nh1.E();

            // tracks loop
            for (Int_t itr2 = 0; itr2 < nTracks; itr2++) {
                track2 = (Track*)branchTrack->At(itr2);
                TLorentzVector cp2;
                cp2.SetPtEtaPhiE(track2->PT, track2->Eta, track2->Phi, track2->P);
                Float_t E2 = cp2.E();
                Float_t cosine_theta = calCosinceTheta(nh1, cp2);
                FWMs_ = updateFWMs(E1, E2, cosine_theta, FWMs_);
            }
            // neutral hadron loop
            for (Int_t inh2 = 0; inh2 < nNeutral; inh2++) {
                neutral2 = (Tower*)branchEFlowNeutralHadron->At(inh2);
                Float_t pt2 = neutral2->ET * pow(neutral2->E * neutral2->E - mK0 * mK0, 0.5) / neutral2->E;
                TLorentzVector nh2;
                nh2.SetPtEtaPhiE(pt2, neutral2->Eta, neutral2->Phi, neutral2->E);
                Float_t E2 = nh2.E();
                Float_t cosine_theta = calCosinceTheta(nh1, nh2);
                FWMs_ = updateFWMs(E1, E2, cosine_theta, FWMs_);
            }
            // photon loop
            for (Int_t iph2 = 0; iph2 < nPhotons; iph2++) {
                photon2 = (Tower*)branchEFlowPhoton->At(iph2);
                TLorentzVector ph2;
                ph2.SetPtEtaPhiE(photon2->ET, photon2->Eta, photon2->Phi, photon2->E);
                Float_t E2 = ph2.E();
                Float_t cosine_theta = calCosinceTheta(nh1, ph2);
                FWMs_ = updateFWMs(E1, E2, cosine_theta, FWMs_);
            }
        }
        // 1st photon loop
        for (Int_t iph1 = 0; iph1 < nPhotons; iph1++) {
            photon1 = (Tower*)branchEFlowPhoton->At(iph1);
            TLorentzVector ph1;
            ph1.SetPtEtaPhiE(photon1->ET, photon1->Eta, photon1->Phi, photon1->E);
            Float_t E1 = ph1.E();

            // tracks loop
            for (Int_t itr2 = 0; itr2 < nTracks; itr2++) {
                track2 = (Track*)branchTrack->At(itr2);
                TLorentzVector cp2;
                cp2.SetPtEtaPhiE(track2->PT, track2->Eta, track2->Phi, track2->P);
                Float_t E2 = cp2.E();
                Float_t cosine_theta = calCosinceTheta(ph1, cp2);
                FWMs_ = updateFWMs(E1, E2, cosine_theta, FWMs_);
            }
            // neutral hadron loop
            for (Int_t inh2 = 0; inh2 < nNeutral; inh2++) {
                neutral2 = (Tower*)branchEFlowNeutralHadron->At(inh2);
                Float_t pt2 = neutral2->ET * pow(neutral2->E * neutral2->E - mK0 * mK0, 0.5) / neutral2->E;
                TLorentzVector nh2;
                nh2.SetPtEtaPhiE(pt2, neutral2->Eta, neutral2->Phi, neutral2->E);
                Float_t E2 = nh2.E();
                Float_t cosine_theta = calCosinceTheta(ph1, nh2);
                FWMs_ = updateFWMs(E1, E2, cosine_theta, FWMs_);
            }
            // photon loop
            for (Int_t iph2 = 0; iph2 < nPhotons; iph2++) {
                photon2 = (Tower*)branchEFlowPhoton->At(iph2);
                TLorentzVector ph2;
                ph2.SetPtEtaPhiE(photon2->ET, photon2->Eta, photon2->Phi, photon2->E);
                Float_t E2 = ph2.E();
                Float_t cosine_theta = calCosinceTheta(ph1, ph2);
                FWMs_ = updateFWMs(E1, E2, cosine_theta, FWMs_);
            }
        }
        /*
        cout << FWMs_.H_EE0 << "; " << FWMs_.H_EE1 << "; " << FWMs_.H_EE2 << "; " << FWMs_.H_EE3 << "; " << FWMs_.H_EE4 << "; "
             << FWMs_.H_EE5 << "; " << FWMs_.H_EE6 << "; " << FWMs_.H_EE7 << "; " << FWMs_.H_EE8 << "; " << FWMs_.H_EE9 << "; " << FWMs_.H_EE10 << endl;
        */
        FWMs->iEvt = i_en;
        FWMs->H_EE0 = FWMs_.H_EE0 / pow(91.2, 2);
        FWMs->H_EE1 = FWMs_.H_EE1 / pow(91.2, 2);
        FWMs->H_EE2 = FWMs_.H_EE2 / pow(91.2, 2);
        FWMs->H_EE3 = FWMs_.H_EE3 / pow(91.2, 2);
        FWMs->H_EE4 = FWMs_.H_EE4 / pow(91.2, 2);
        FWMs->H_EE5 = FWMs_.H_EE5 / pow(91.2, 2);
        FWMs->H_EE6 = FWMs_.H_EE6 / pow(91.2, 2);
        FWMs->H_EE7 = FWMs_.H_EE7 / pow(91.2, 2);
        FWMs->H_EE8 = FWMs_.H_EE8 / pow(91.2, 2);
        FWMs->H_EE9 = FWMs_.H_EE9 / pow(91.2, 2);
        FWMs->H_EE10 = FWMs_.H_EE10 / pow(91.2, 2);
        tr.Fill();
    }
    tr.Write();
    fea.Close();
}
