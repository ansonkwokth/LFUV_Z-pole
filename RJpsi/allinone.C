#include <cmath>
#include <iostream>
#include <random>
#include <vector>
using namespace std;

#include "FeaturesClass.C"
#include "FinalStatesClass.C"
#include "Geometry.C"
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

/*
//Load Delphes for reading the data
#ifdef __CLING__
R__ADD_LIBRARY_PATH($DELPHES)
R__LOAD_LIBRARY(libDelphes)

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#else
class ExRootResult;
class ExRootTreeReader;
#endif

*/

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

// Find final states index {{{
iFinalStates FindFinalStatesIndex(TClonesArray* branchTrack, Int_t nMu = 3) {
    const Float_t mMu = 0.1057;
    iFinalStates iFinalStatesIndexes;
    Int_t nTracks = branchTrack->GetEntries();
    Track* track1;
    Track* track2;

    Int_t numMuTot = 0;
    for (int it = 0; it < nTracks; it++) {
        track1 = (Track*)branchTrack->At(it);
        if (abs(track1->PID) == 13) {
            numMuTot += 1;
        }
    }
    if (numMuTot != nMu) {
        return iFinalStatesIndexes;
    }

    // finding positive muon decay vertex
    Int_t vert[nTracks];
    Int_t iLoc = 0;
    for (int it = 0; it < nTracks; it++) {
        track1 = (Track*)branchTrack->At(it);
        if (track1->PID == -13 && Length(track1->X, track1->Y, track1->Z) > 0) {
            vert[iLoc] = it;
            iLoc += 1;
        }
    }
    if (iLoc < 1) {
        return iFinalStatesIndexes;
    }

    Int_t foundMuNeg = 0;
    Int_t iMuPos, iMuNeg;
    Float_t diffM = 99999;
    // finding same vertex oppositive charge pion
    for (int iar = 0; iar < iLoc; iar++) {
        track1 = (Track*)branchTrack->At(vert[iar]);
        // cout<< iar  <<endl;
        // Int_t foundMuNegi = 0;
        for (int it2 = 0; it2 < nTracks; it2++) {
            // same vertex
            track2 = (Track*)branchTrack->At(it2);
            if (track2->PID == 13 && track1->X == track2->X && track1->Y == track2->Y && track1->Z == track2->Z && track1->Charge + track2->Charge == 0) {
                TLorentzVector mu1, mu2, Hc;

                mu1.SetPtEtaPhiM(track1->PT, track1->Eta, track1->Phi, mMu);
                mu2.SetPtEtaPhiM(track2->PT, track2->Eta, track2->Phi, mMu);
                Hc = mu1 + mu2;
                if (abs(Hc.M() - 3.09) < diffM) {
                    // cout<< Hc.M() <<endl;
                    diffM = abs(Hc.M() - 3.09);
                    iMuPos = vert[iar];
                    iMuNeg = it2;
                    // foundMuNegi = 1;
                    foundMuNeg = 1;
                } else {
                    continue;
                }
                /*
                iMuPos = vert[iar];
                iMuNeg = it2;
                foundMuNeg = 1;
                break;
                */
            }
        }
        // if (foundMuNeg == 1) { break; }
    }

    if (not(foundMuNeg == 1)) {
        return iFinalStatesIndexes;
    }
    iFinalStatesIndexes.iMuPos = iMuPos;
    iFinalStatesIndexes.iMuNeg = iMuNeg;
    iFinalStatesIndexes.foundFromC = 1;
    return iFinalStatesIndexes;
}

//}}}

//---------- {{{
Int_t findMuIndex(iFinalStates iFS, TClonesArray* branchTrack) {
    Int_t nTracks = branchTrack->GetEntries();
    Track* trackMu;
    // find muon
    Int_t foundMu = 0;
    Int_t iMu;
    for (int it = 0; it < nTracks; it++) {
        trackMu = (Track*)branchTrack->At(it);
        if (abs(trackMu->PID) == 13 && it != iFS.iMuPos && it != iFS.iMuNeg) {
            iMu = it;
            foundMu += 1;
        }
    }
    if (foundMu != 1) {
        return 99999;
    }
    return iMu;
}
//}}}

//----------i {{{
TLorentzVector
reconstructHc4Momentum(TVector3 v3B, iFinalStates iFS, TClonesArray* branchTrack, TClonesArray* branchEFlowPhoton, TClonesArray* branchEFlowNeutralHadron) {
    const Float_t mK0 = 0.497661;
    const Float_t mBc = 6.2749;
    TLorentzVector signalHemi, recoilHemi, pRem;

    Int_t nTracks = branchTrack->GetEntries();
    Int_t nEFlowPhotons = branchEFlowPhoton->GetEntries();
    Int_t nEFlowNeutralHadrons = branchEFlowNeutralHadron->GetEntries();

    Track* track;
    Tower* eflowph;
    Tower* eflownh;

    // charged particles
    if (nTracks > 0) {
        for (int it = 0; it < nTracks; it++) {
            track = (Track*)branchTrack->At(it);
            TLorentzVector chargedPi;
            chargedPi.SetPtEtaPhiE(track->PT, track->Eta, track->Phi, track->P);
            if (chargedPi.Px() * v3B.X() + chargedPi.Py() * v3B.Y() + chargedPi.Pz() * v3B.Z() > 0) {
                signalHemi += chargedPi;
                if (it != iFS.iMuPos && it != iFS.iMuNeg && it != iFS.iMu) {
                    pRem += chargedPi;
                }
            } else {
                recoilHemi += chargedPi;
            }
        }
    }

    // photon
    if (nEFlowPhotons > 0) {
        for (int iefp = 0; iefp < nEFlowPhotons; iefp++) {
            eflowph = (Tower*)branchEFlowPhoton->At(iefp);
            TLorentzVector eflowPhoi;
            eflowPhoi.SetPtEtaPhiE(eflowph->ET, eflowph->Eta, eflowph->Phi, eflowph->E);
            if (eflowPhoi.Px() * v3B.X() + eflowPhoi.Py() * v3B.Y() + eflowPhoi.Pz() * v3B.Z() > 0) {
                signalHemi += eflowPhoi;
                pRem += eflowPhoi;
            } else {
                recoilHemi += eflowPhoi;
            }
        }
    }

    // neutral hadrons
    if (nEFlowNeutralHadrons > 0) {
        for (int iefnh = 0; iefnh < nEFlowNeutralHadrons; iefnh++) {
            eflownh = (Tower*)branchEFlowNeutralHadron->At(iefnh);
            Float_t pt = eflownh->ET * pow(eflownh->E * eflownh->E - mK0 * mK0, 0.5) / eflownh->E;
            TLorentzVector eflowNeuHi;
            eflowNeuHi.SetPtEtaPhiE(pt, eflownh->Eta, eflownh->Phi, eflownh->E);
            if (eflowNeuHi.Px() * v3B.X() + eflowNeuHi.Py() * v3B.Y() + eflowNeuHi.Pz() * v3B.Z() > 0) {
                signalHemi += eflowNeuHi;
                pRem += eflowNeuHi;
            } else {
                recoilHemi += eflowNeuHi;
            }
        }
    }

    Float_t EB = (91.0 * 91.0 + signalHemi.M() * signalHemi.M() - recoilHemi.M() * recoilHemi.M()) / (2 * 91.0) - pRem.E();

    Float_t pB = pow(EB * EB - mBc * mBc, 0.5);
    TLorentzVector B;
    Float_t ratio = pB / (pow(v3B.X() * v3B.X() + v3B.Y() * v3B.Y() + v3B.Z() * v3B.Z(), 0.5));
    B.SetPxPyPzE(v3B.X() * ratio, v3B.Y() * ratio, v3B.Z() * ratio, EB);

    return B;
}
//}}}

//----------{{{
vector<Int_t>
findMisIDPi(iFinalStates iFS, TClonesArray* branchTrack) {
    vector<Int_t> iPiS;
    Int_t nTracks = branchTrack->GetEntries();
    Track* track;
    for (int it = 0; it < nTracks; it++) {
        if (it == iFS.iMuPos || it == iFS.iMuNeg || it == iFS.iMu) {
            continue;
        }
        track = (Track*)branchTrack->At(it);
        if (abs(track->PID) != 211) {
            continue;
        }
        if (track->X == 0 && track->Y == 0 && track->Z == 0) {
            continue;
        }
        iPiS.push_back(it);
    }
    return iPiS;
}
//}}}

//----------{{{
Float_t
cal_pPerp(TVector3 v3B, TLorentzVector Hc, TLorentzVector Mu) {
    TLorentzVector HcMu;
    HcMu = Hc + Mu;
    Float_t proj = (HcMu.Px() * v3B.X() + HcMu.Py() * v3B.Y() + HcMu.Pz() * v3B.Z()) / pow(Length(v3B.X(), v3B.Y(), v3B.Z()), 2);
    Float_t pPerp = Length(HcMu.Px() - proj * v3B.X(), HcMu.Py() - proj * v3B.Y(), HcMu.Pz() - proj * v3B.Z());
    return pPerp;
}
///}}}

//-----{{{
Float_t cal_mCorr(TVector3 v3B, TLorentzVector Hc, TLorentzVector Mu) {
    Float_t pPerp = cal_pPerp(v3B, Hc, Mu);
    TLorentzVector HcMu;
    HcMu = Hc + Mu;
    Float_t mCorr = pow(HcMu.M() * HcMu.M() + pPerp * pPerp, 0.5) + pPerp;
    return mCorr;
}
//}}}

//----------{{{
struct impactParams {
    Float_t D0Max = 0;
    Float_t D0Sum = 0;
    Float_t DzMax = 0;
    Float_t DzSum = 0;
};
//}}}

//----------{{{
impactParams FindImpactParams(TClonesArray* branchTrack, iFinalStates iFS, TVector3 v3B) {
    impactParams impParams;

    Int_t nTracks = branchTrack->GetEntries();
    Track* track;
    for (int it = 0; it < nTracks; it++) {
        track = (Track*)branchTrack->At(it);
        if (it == iFS.iMuPos || it == iFS.iMuNeg || it == iFS.iMu) {
            continue;
        }
        TLorentzVector impTr;
        impTr.SetPtEtaPhiE(track->PT, track->Eta, track->Phi, track->P);
        if ((impTr.Px() * v3B.X() + impTr.Py() * v3B.Y() + impTr.Pz() * v3B.Z() > 0) && (track->X * v3B.X() + track->Y * v3B.Y() + track->Z * v3B.Z() > 0)) {
            impParams.D0Sum += abs(track->D0);
            impParams.DzSum += abs(track->DZ);
            if (abs(track->D0) > impParams.D0Max) {
                impParams.D0Max = abs(track->D0);
            }
            if (abs(track->DZ) > impParams.DzMax) {
                impParams.DzMax = abs(track->DZ);
            }
        }
    }
    return impParams;
}

//}}}

//----------{{{
struct isolationVars {
    Float_t ENeutral03 = 0;
    Float_t ENeutral06 = 0;
    Float_t ENeutral03Hadron = 0;
    Float_t ENeutral06Hadron = 0;
    Float_t ENeutral03Photon = 0;
    Float_t ENeutral06Photon = 0;

    Float_t ECharge03 = 0;
    Float_t ECharge06 = 0;
    Float_t ECharge03PV = 0;
    Float_t ECharge06PV = 0;
    Float_t ECharge03DV = 0;
    Float_t ECharge06DV = 0;
};
//}}}

//---------- {{{
isolationVars
FindIsolationVars(TClonesArray* branchTrack, TClonesArray* branchEFlowNeutralHadron, TClonesArray* branchEFlowPhoton, iFinalStates iFS, TLorentzVector Hc) {
    isolationVars isoVars;

    const Float_t mK0 = 0.497661;
    TLorentzVector B03Neu, B06Neu, B03NeuHad, B06NeuHad, B03NeuPho, B06NeuPho;
    TLorentzVector B03Chg, B06Chg, B03ChgPV, B06ChgPV, B03ChgDV, B06ChgDV;

    Int_t nEFlowNeutralHadrons = branchEFlowNeutralHadron->GetEntries();
    Int_t nEFlowPhotons = branchEFlowPhoton->GetEntries();
    Int_t nTracks = branchTrack->GetEntries();

    Tower* eflowNeuHad;
    Tower* eflowPho;
    Track* track;

    for (int in = 0; in < nEFlowNeutralHadrons; in++) {
        eflowNeuHad = (Tower*)branchEFlowNeutralHadron->At(in);
        Float_t pt = eflowNeuHad->ET * pow(eflowNeuHad->E * eflowNeuHad->E - mK0 * mK0, 0.5) / eflowNeuHad->E;
        TLorentzVector NeuHadi;
        NeuHadi.SetPtEtaPhiE(pt, eflowNeuHad->Eta, eflowNeuHad->Phi, eflowNeuHad->E);
        if (Angle(NeuHadi.Px(), NeuHadi.Py(), NeuHadi.Pz(), Hc.Px(), Hc.Py(), Hc.Pz()) < 0.3) {
            B03Neu += NeuHadi;
            B03NeuHad += NeuHadi;
        }
        if (Angle(NeuHadi.Px(), NeuHadi.Py(), NeuHadi.Pz(), Hc.Px(), Hc.Py(), Hc.Pz()) < 0.6) {
            B06Neu += NeuHadi;
            B06NeuHad += NeuHadi;
        }
    }

    for (int iph = 0; iph < nEFlowPhotons; iph++) {
        eflowPho = (Tower*)branchEFlowPhoton->At(iph);
        TLorentzVector Phoi;
        Phoi.SetPtEtaPhiE(eflowPho->ET, eflowPho->Eta, eflowPho->Phi, eflowPho->E);
        if (Angle(Phoi.Px(), Phoi.Py(), Phoi.Pz(), Hc.Px(), Hc.Py(), Hc.Pz()) < 0.3) {
            B03Neu += Phoi;
            B03NeuPho += Phoi;
        }
        if (Angle(Phoi.Px(), Phoi.Py(), Phoi.Pz(), Hc.Px(), Hc.Py(), Hc.Pz()) < 0.6) {
            B06Neu += Phoi;
            B06NeuPho += Phoi;
        }
    }
    for (int it = 0; it < nTracks; it++) {
        if (it == iFS.iMuPos || it == iFS.iMuNeg || it == iFS.iMu) {
            continue;
        }
        track = (Track*)branchTrack->At(it);
        TLorentzVector Tri;
        Tri.SetPtEtaPhiE(track->PT, track->Eta, track->Phi, track->P);
        if (Angle(Tri.Px(), Tri.Py(), Tri.Pz(), Hc.Px(), Hc.Py(), Hc.Pz()) < 0.3) {
            B03Chg += Tri;
            if (track->X == 0 && track->Y == 0 && track->Z == 0) {
                B03ChgPV += Tri;
            } else {
                B03ChgDV += Tri;
            }
        }
        if (Angle(Tri.Px(), Tri.Py(), Tri.Pz(), Hc.Px(), Hc.Py(), Hc.Pz()) < 0.6) {
            B06Chg += Tri;
            if (track->X == 0 && track->Y == 0 && track->Z == 0) {
                B06ChgPV += Tri;
            } else {
                B06ChgDV += Tri;
            }
        }
    }

    isoVars.ENeutral03 = B03Neu.E();
    isoVars.ENeutral06 = B06Neu.E();
    isoVars.ENeutral03Hadron = B03NeuHad.E();
    isoVars.ENeutral06Hadron = B06NeuHad.E();
    isoVars.ENeutral03Photon = B03NeuPho.E();
    isoVars.ENeutral06Photon = B06NeuPho.E();

    isoVars.ECharge03 = B03Chg.E();
    isoVars.ECharge06 = B06Chg.E();
    isoVars.ECharge03PV = B03ChgPV.E();
    isoVars.ECharge06PV = B06ChgPV.E();
    isoVars.ECharge03DV = B03ChgDV.E();
    isoVars.ECharge06DV = B06ChgDV.E();
    // cout << isoVars.ECharge06DV<<endl;
    return isoVars;
}
//}}}

//{{{
TLorentzVector
reconstructK0S(TClonesArray* branchTrack, iFinalStates iFS, double noise) {
    // reconstruct K0S -> pi pi
    TLorentzVector K0S;
    K0S.SetPtEtaPhiE(99999, 99999, 99999, 99999);

    // inject noise, in track location
    std::default_random_engine genertator;
    std::normal_distribution<double> distribution(0, noise / pow(3, 0.5));

    // store the pions cadidates
    vector<Int_t> iPosS;
    vector<Int_t> iNegS;
    Int_t nTracks = branchTrack->GetEntries();
    Track* track;
    Track* trackPos;
    Track* trackNeg;
    for (int it = 0; it < nTracks; it++) {
        if (it == iFS.iMuPos || it == iFS.iMuNeg || it == iFS.iMu) {
            continue;
        }
        track = (Track*)branchTrack->At(it);
        if (track->X == 0 && track->Y == 0 && track->Z == 0) {
            continue;
        }
        if (track->PID == 211) {
            iPosS.push_back(it);
        }
        if (track->PID == -211) {
            iNegS.push_back(it);
        }
    }

    Float_t mK0SDelta = 99999;
    TVector3 v3K0S;
    Int_t savediPos, savediNeg;
    for (int ipos : iPosS) {
        trackPos = (Track*)branchTrack->At(ipos);
        for (int ineg : iNegS) {
            trackNeg = (Track*)branchTrack->At(ineg);
            if (not(trackPos->X == trackNeg->X && trackPos->Y && trackNeg->Y && trackPos->Z == trackNeg->Z)) {
                continue;
            }

            Float_t XK0S = trackPos->X + distribution(genertator);
            Float_t YK0S = trackPos->Y + distribution(genertator);
            Float_t ZK0S = trackPos->Z + distribution(genertator);

            if (Length(XK0S, YK0S, ZK0S) <= 10) {
                continue;
            }
            // cout << "Pos pi vertex length: " << Length(trackPos->X,
            // trackPos->Y, trackPos->Z) << endl; cout << "Neg pi vertex length:
            // " << Length(trackNeg->X, trackNeg->Y, trackNeg->Z) << endl;
            TLorentzVector piPos;
            piPos.SetPtEtaPhiM(trackPos->PT, trackPos->Eta, trackPos->Phi, 0.13957);
            TLorentzVector piNeg;
            piNeg.SetPtEtaPhiM(trackNeg->PT, trackNeg->Eta, trackNeg->Phi, 0.13957);
            TLorentzVector K0Si;
            K0Si = piPos + piNeg;
            if (abs(K0Si.M() - 0.49761) < mK0SDelta) {
                savediPos = ipos;
                savediNeg = ineg;
                v3K0S.SetXYZ(XK0S, YK0S, ZK0S);
                K0S = K0Si;
                mK0SDelta = abs(K0Si.M() - 0.49761);
            }
            // cout << abs(K0Si.M() - 0.49761) << endl;
        }
    }
    // cout << abs(K0S.M() - 0.49761) << endl;
    // cout << Length(v3K0S.X(), v3K0S.Y(), v3K0S.Z()) << endl;
    // cout <<endl;

    Float_t disTargetTr = 99999;  // closest distance between a track and the target
    Track* trackOther;
    TLorentzVector target = K0S;
    TVector3 v3Target = v3K0S;

    for (int it = 0; it < nTracks; it++) {
        if (it == iFS.iMuPos || it == iFS.iMuNeg || it == iFS.iMu || it == savediNeg || it == savediPos) {
            continue;
        }
        trackOther = (Track*)branchTrack->At(it);
        if (trackOther->X == 0 && trackOther->Y == 0 && trackOther->Z == 0) {
            continue;
        }

        TLorentzVector otherTrack;
        otherTrack.SetPtEtaPhiE(trackOther->PT, trackOther->Eta, trackOther->Phi, trackOther->P);
        if (otherTrack.Px() * target.Px() + otherTrack.Py() * target.Py() + otherTrack.Pz() * target.Pz() <= 0) {
            continue;
        }

        Float_t XTr = trackOther->X + distribution(genertator);
        Float_t YTr = trackOther->Y + distribution(genertator);
        Float_t ZTr = trackOther->Z + distribution(genertator);

        Float_t L, s1, s2;
        distance_2lines(XTr, YTr, ZTr, otherTrack.Px(), otherTrack.Py(), otherTrack.Pz(), v3Target.X(), v3Target.Y(), v3Target.Z(), target.Px(), target.Py(), target.Pz(), &L, &s1, &s2);
        if (L < disTargetTr) {
            disTargetTr = L;
        }
    }
    // veto K0S
    Float_t disK0STr = disTargetTr;
    // disK0STr = closestTrack(iFS, K0S, v3K0S, noise, branchTrack);
    // cout << "disTargetTr " << disK0STr << endl;
    if (disK0STr <= 0.02) {
        K0S.SetPtEtaPhiE(99999, 99999, 99999, 99999);
    }

    return K0S;
}
//}}}

void allinone(const string type, const Float_t noise_ = 10, const Bool_t save = true, Int_t num_test = 0) {
    cout << "\n\n\n\n\n\n\n\n\n\n"
         << endl;
    string typeName;
    const char* inputFile;
    const char* outputFile;
    const Float_t noise = noise_ * 0.001;
    if (type == "s1") {
        typeName = "Bc->Jpsi tau nu. ";
        cout << "Bc->Jpsi tau nu. " << endl;
        inputFile = "./Bcjpsitaunu_50m.root";
        if (noise_ == 10) {
            outputFile = "./features/JpsiTauNu_10Noise_NoVeto.root";
        } else if (noise_ == 20) {
            outputFile = "./features/JpsiTauNu_20Noise_NoVeto.root";
        }
    } else if (type == "s2") {
        cout << "Bc->Jpsi mu nu. " << endl;
        inputFile = "./BcJpsimunu0-2.root";
        if (noise_ == 10) {
            outputFile = "./features/JpsiMuNu_10Noise_NoVeto.root";
        } else if (noise_ == 20) {
            outputFile = "./features/JpsiMuNu_20Noise_NoVeto.root";
        }
    } else if (type == "b1") {
        cout << "Comb+Cascade Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        if (noise_ == 10) {
            outputFile = "./features/RJpsiCombCascade_10Noise_NoVeto.root";
        } else if (noise_ == 20) {
            outputFile = "./features/RJpsiCombCascade_20Noise_NoVeto.root";
        }
    } else if (type == "b2") {
        cout << "Comb Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        if (noise_ == 10) {
            outputFile = "./features/RJpsiComb_10Noise_NoVeto.root";
        } else if (noise_ == 20) {
            outputFile = "./features/RJpsiComb_20Noise_NoVeto.root";
        }
    } else if (type == "b3") {
        cout << "Cascade Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        if (noise_ == 10) {
            outputFile = "./features/RJpsiCascade_10Noise_NoVeto.root";
        } else if (noise_ == 20) {
            outputFile = "./features/RJpsiCascade_20Noise_NoVeto.root";
        }
    } else if (type == "b4") {
        cout << "Inclusive Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        if (noise_ == 10) {
            outputFile = "./features/RJpsiInclusive_10Noise_NoVeto.root";
        } else if (noise_ == 20) {
            outputFile = "./features/RJpsiInclusive_20Noise_NoVeto.root";
        }
    } else if (type == "b5") {
        cout << "MisID Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        if (noise_ == 10) {
            outputFile = "./features/RJpsiMisID_10Noise_NoVeto.root";
        } else if (noise_ == 20) {
            outputFile = "./features/RJpsiMisID_20Noise_NoVeto.root";
        }
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

    // inject noise
    std::default_random_engine genertator;
    std::normal_distribution<double> distribution(0, noise / pow(3, 0.5));

    // defining constants
    const Float_t mK = 0.493677;
    const Float_t mK0 = 0.497661;
    const Float_t mP = 0.938272;
    const Float_t mPi = 0.13957;
    const Float_t mMu = 0.1057;
    const Float_t mBc = 6.2749;
    const Float_t mTau = 1.777;
    const Float_t mLambdab = 5.61960;

    GenParticle* particle;
    Track* track;
    Tower* eflowphoton;
    Tower* eflowneutralhadron;

    // Count number of targeted events
    Int_t nEvt = 0;
    Int_t nFoundTracks = 0;
    Int_t nMu = 0;
    Int_t nSameDir = 0;
    Int_t nHcMass = 0;
    Int_t nMuPt = 0;
    Int_t nVert = 0;
    Int_t num = 0;

    // output daata
    if (not save) {
        outputFile = "dummy.root";
    }
    TFile fea(outputFile, "recreate");
    TTree tr("t", "features");
    Features* features = new Features;
    // Features features;
    tr.Branch("features", &features);

    cout << "\nReconstructing: " << typeName << endl;
    cout << "Reading from:\t" << inputFile << endl;
    cout << "Injected noise:\t" << noise_ << " microns (" << noise << "mm)" << endl;
    cout << "Save: " << save << endl;
    cout << "Writing to:\t" << outputFile << endl;
    cout << "# of events:\t" << numberOfEntries << "\n"
         << endl;
    if (num_test == 0) {
        num_test = numberOfEntries;
    }

    for (Int_t i_en = 0; i_en < numberOfEntries; i_en++) {
        // for (Int_t i_en=0; i_en<num_test; i_en++) {
        // if (i_en % 1000 == 0) {cout << "\tReconstruction Progress: " << i_en
        // << "/" << numberOfEntries; } if ((i_en % 1000) == 0) {cout <<
        // "\rReconstruction Progress: " << i_en << "/" << numberOfEntries; }
        if ((i_en % 100000) == 0) {
            cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << endl;
        }
        if (i_en >= num_test) {
            break;
        }
        treeReader->ReadEntry(i_en);  // reading the entry

        //==========================================================================
        //===============   Classifying in event type in truth level
        //===============
        //==========================================================================
        iFinalStates iFSTrue;
        TLorentzVector BTrue, HcTrue, muTrue;
        TVector3 v3HcTrue, v3MuTrue;
        Int_t passing = 0;

        // cout<<"----------------Event: " << i_en <<endl;
        // cout<<  "" <<endl;
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

        if (passing == 0) {
            continue;
        }

        nEvt += 1;
        // cout<< "...............Identified: "<<i_en<<"; " <<nEvt<<endl;
        //==========================================================================
        //===============   Finding the correspdoning final states
        //===============
        //==========================================================================
        Int_t nMu = 3;
        if (type == "b5") {
            nMu = 2;
        }
        iFinalStates iFS = FindFinalStatesIndex(branchTrack, nMu = nMu);
        if (iFS.foundFromC == 0) {
            continue;
        }
        nFoundTracks += 1;

        vector<Int_t> muCandidates;
        if (type == "b5") {
            muCandidates = findMisIDPi(iFS, branchTrack);
        } else {
            muCandidates.push_back(findMuIndex(iFS, branchTrack));
        }

        Int_t nRecoMu = 0;
        for (Int_t imu : muCandidates) {
            // cout<< imu <<endl;
            if (imu == 99999) {
                continue;
            }
            // cout << iloop << endl;
            iFS.iMu = imu;

            //==========================================================================
            //==================   Defining 4 Momentum & 3 Positions  ==================
            //==========================================================================
            Track* muPosTrack = (Track*)branchTrack->At(iFS.iMuPos);
            Track* muNegTrack = (Track*)branchTrack->At(iFS.iMuNeg);
            Track* muTrack = (Track*)branchTrack->At(iFS.iMu);

            TLorentzVector muPos;
            muPos.SetPtEtaPhiM(muPosTrack->PT, muPosTrack->Eta, muPosTrack->Phi, mMu);
            TLorentzVector muNeg;
            muNeg.SetPtEtaPhiM(muNegTrack->PT, muNegTrack->Eta, muNegTrack->Phi, mMu);
            TLorentzVector mu;
            mu.SetPtEtaPhiM(muTrack->PT, muTrack->Eta, muTrack->Phi, mMu);

            if (muPos.Px() * muNeg.Px() + muPos.Py() * muNeg.Py() + muPos.Pz() * muNeg.Pz() <= 0) {
                continue;
            }
            if (muPos.Px() * mu.Px() + muPos.Py() * mu.Py() + muPos.Pz() * mu.Pz() <= 0) {
                continue;
            }
            if (muNeg.Px() * mu.Px() + muNeg.Py() * mu.Py() + muNeg.Pz() * mu.Pz() <= 0) {
                continue;
            }

            nSameDir += 1;

            TLorentzVector Hc;
            Hc = muPos + muNeg;

            Float_t dxC = distribution(genertator);
            Float_t dyC = distribution(genertator);
            Float_t dzC = distribution(genertator);
            Float_t dxMu = distribution(genertator);
            Float_t dyMu = distribution(genertator);
            Float_t dzMu = distribution(genertator);
            TVector3 v3CNoise(dxC, dyC, dzC);
            TVector3 v3MuNoise(dxMu, dyMu, dzMu);
            TVector3 v3C(muPosTrack->X, muPosTrack->Y, muPosTrack->Z);
            TVector3 v3Mu(muTrack->X, muTrack->Y, muTrack->Z);
            v3Mu += v3MuNoise;
            v3C += v3CNoise;

            // cout<<  muTrack->X<<"/"<<muPosTrack->X << "; " << muTrue.X()
            // <<endl;

            //==========================================================================
            //========================   Vetoing muon and Hc
            //=========================
            //==========================================================================

            // veto muon
            Float_t disMuTr;
            disMuTr = closestTrack(iFS, mu, v3Mu, noise, branchTrack);
            // disMuTr = closestTrack_debug(iFS, mu, v3Mu, noise, branchTrack,
            // branchParticle); cout<<  "" <<endl; cout <<
            // "\n===========================================Event: " << i_en <<
            // endl; cout << " ..." << iFS.iMuPos << "; " << iFS.iMuNeg << "; "
            // << iFS.iMu << endl; cout << "Closest Distance: " << disMuTr
            // <<endl;
            // if (disMuTr < 0.02) {
            //    continue;
            //}
            nMu += 1;

            // veto Hc
            Float_t disHcTr;
            disHcTr = closestTrack(iFS, Hc, v3C, noise, branchTrack);
            // if (disHcTr < 0.02) { continue; }

            //==========================================================================
            //=============================   Apply cuts   =============================
            //==========================================================================
            // cout<<  "" <<endl;
            // cout<<  "=====" << i_en << "; " << nVert<<endl;
            if (not(abs(Hc.M() - 3.09) < 0.055 / 2 && Hc.Pt() > 2.0 / 2 && muPos.P() > 5.0 / 2 && muNeg.P() > 5.0 / 2 && (muPos.Pt() > 1.5 / 2 || muNeg.Pt() > 1.5 / 2))) {
                continue;
            }
            nHcMass += 1;

            if (not(mu.Pt() > 0.75 / 2 && mu.P() > 3.0 / 2 && mu.P() < 100 * 2 && (mu.Px() * Hc.Px() + mu.Py() * Hc.Py() + mu.Pz() * Hc.Pz() > 0))) {
                continue;
            }

            nMuPt += 1;

            //==========================================================================
            //=======================   Deduce B decay vertex
            //========================
            //==========================================================================
            if (Length(v3C.X(), v3C.Y(), v3C.Z()) < 0.1) {
                continue;
            }
            nVert += 1;
            // cout<<  Length(v3C.X(), v3C.Y(), v3C.Z())<<endl;
            // deduced location
            Float_t X = v3C.X();
            Float_t Y = v3C.Y();
            Float_t Z = v3C.Z();
            TVector3 v3B(X, Y, Z);

            /*
            if (Length(v3C.X(), v3C.Y(), v3C.Z()) < 0.05) {
                continue;
            }
            nVert += 1;
            */
            TLorentzVector B;
            B = reconstructHc4Momentum(v3B, iFS, branchTrack, branchEFlowPhoton, branchEFlowNeutralHadron);
            if (isnan(B.P())) {
                continue;
            }
            if (B.Px() == 99999 && B.Py() == 99999 && B.Pz() == 99999 && B.E() == 99999) {
                continue;
            }

            TLorentzVector HcMu;
            HcMu = muPos + muNeg + mu;
            if (HcMu.M() > mBc) {
                continue;
            }

            //==========================================================================
            //==============================   Feature    ==============================
            //==========================================================================
            features->iEvt = i_en;
            // q2
            TLorentzVector q;
            q = B - Hc;
            Float_t q2 = q.E() * q.E() - q.P() * q.P();
            if (q2 < -10 || q2 > 15) {
                continue;
            }
            features->q2 = q2;
            // miss2
            TLorentzVector miss;
            miss = B - Hc - mu;
            features->miss2 = miss.E() * miss.E() - miss.P() * miss.P();
            // pB, EB
            features->pB = B.P();
            features->EB = B.E();
            // pHc, EHc
            features->pHc = Hc.P();
            features->EHc = Hc.E();
            // pMu, EMu
            features->pMu = mu.P();
            features->EMu = mu.E();
            // distances
            features->sMinMuHcVert = distance_linepoint(v3C, v3Mu, mu);
            features->sMinMuTr = disMuTr;
            features->sMinHcTr = disHcTr;
            features->sPVHc = Length(v3C.X(), v3C.Y(), v3C.Z());
            features->mHcMu = HcMu.M();
            // p_perp
            features->pPerp = cal_pPerp(v3B, Hc, mu);
            // m_corr
            features->mCorr = cal_mCorr(v3B, Hc, mu);

            // impact parameters
            impactParams impParams;
            impParams = FindImpactParams(branchTrack, iFS, v3B);
            features->D0Max = impParams.D0Max;
            features->D0Sum = impParams.D0Sum;
            features->DzMax = impParams.DzMax;
            features->DzSum = impParams.DzSum;

            // isolation variables
            isolationVars isoVars;
            isoVars = FindIsolationVars(branchTrack, branchEFlowNeutralHadron, branchEFlowPhoton, iFS, Hc);

            features->ENeutral03 = isoVars.ENeutral03;
            features->ENeutral06 = isoVars.ENeutral06;
            features->ENeutral03Hadron = isoVars.ENeutral03Hadron;
            features->ENeutral06Hadron = isoVars.ENeutral06Hadron;
            features->ENeutral03Photon = isoVars.ENeutral03Photon;
            features->ENeutral06Photon = isoVars.ENeutral06Photon;

            features->ECharge03 = isoVars.ECharge03;
            features->ECharge06 = isoVars.ECharge06;
            features->ECharge03PV = isoVars.ECharge03PV;
            features->ECharge06PV = isoVars.ECharge06PV;
            features->ECharge03DV = isoVars.ECharge03DV;
            features->ECharge06DV = isoVars.ECharge06DV;

            TLorentzVector qTrue;
            qTrue = BTrue - HcTrue;
            TLorentzVector missTrue;
            missTrue = BTrue - HcTrue - muTrue;
            Float_t q2True = qTrue.E() * qTrue.E() - qTrue.P() * qTrue.P();
            Float_t miss2True = missTrue.E() * missTrue.E() - missTrue.P() * missTrue.P();

            features->q2True = q2True;
            features->miss2True = miss2True;
            features->EBTrue = BTrue.E();
            features->pBTrue = BTrue.P();
            // cout << " ...." << endl;
            // cout << "ien: "<<i_en <<"; "<<  q2 << "; " << q2True<<endl;
            // cout << X << "; " << Y << "; " << Z << endl;
            // cout<< mu.Px() << "; " << mu.Py() << "; " << mu.Pz() <<endl;
            // cout<< muPos.Px() << "; " << muPos.Py() << "; " << muPos.Pz()
            // <<endl; cout<< muNeg.Px() << "; " << muNeg.Py() << "; " <<
            // muNeg.Pz() <<endl; cout<< Hc.Px() << "; " << Hc.Py() << "; " <<
            // Hc.Pz() <<endl; cout<< B.Px() << "; " << B.Py() << "; " << B.Pz()
            // <<endl; cout<< q.Px() << "; " << q.Py() << "; " << q.Pz() <<endl;
            /*
            cout << features->ENeutral03 <<"; " <<
            features->ENeutral06  << "; " <<
            features->ENeutral03Hadron   << "; " <<
            features->ENeutral06Hadron   << "; " <<
            features->ENeutral03Photon   << "; " <<
            features->ENeutral06Photon   << endl;

            cout << features->ECharge03   << "; " <<
            features->ECharge06   << "; " <<
            features->ECharge03PV   << "; " <<
            features->ECharge06PV   << "; " <<
            features->ECharge03DV   << "; " <<
            features->ECharge06DV   << ";" << endl;
            */
            TLorentzVector K0S;
            K0S = reconstructK0S(branchTrack, iFS, noise);
            TLorentzVector K0SHcMu;
            K0SHcMu = K0S + Hc + mu;
            features->mK0SHcMu = K0SHcMu.M();
            features->pK0S = K0S.P();

            Float_t sMinMuHcVertTrue;
            sMinMuHcVertTrue = distance_linepoint(v3HcTrue, v3MuTrue, muTrue);
            features->sMinMuHcVertTrue = sMinMuHcVertTrue;
            tr.Fill();
            // cout << features->sMinMuHcVert << endl;
            //  cout<<features->q2  << "; " << features->q2True<<endl;
            num += 1;
            nRecoMu += 1;
        }
        // cout << nRecoMu << endl;
        //  cout << endl;
    }

    cout << endl;
    cout << endl;
    cout << endl;
    cout << "Number of Events: \t\t" << nEvt << endl;
    cout << "Matched vertex: \t\t" << nFoundTracks << endl;
    cout << "All in same direction: \t\t" << nSameDir << endl;
    cout << "Matched muon: \t\t\t" << nMu << endl;
    cout << "Hc in mass range: \t\t" << nHcMass << endl;
    cout << "Mu PT cut: \t\t\t" << nMuPt << endl;
    cout << "Deduced vertex cut: \t\t" << nVert << endl;
    cout << "Number of Reconstrcutions: \t" << num << " / " << nEvt << endl;

    if (save) {
        tr.Write();
        fea.Close();
        cout << "Writing to:\t" << outputFile << endl;
    }
}

