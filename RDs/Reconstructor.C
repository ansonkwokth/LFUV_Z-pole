#pragma once
#include <iostream>
using namespace std;
#include <vector>

#include "FinalStatesClass.C"
#include "Geometry.C"
#include "Rtypes.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

// tage the final states particles (decay from c-hadron)
iFinalStates FindFinalStatesIndex(TClonesArray* branchTrack) {
    iFinalStates iFinalStatesIndexes;
    Int_t nTracks = branchTrack->GetEntries();
    Track* track1;
    Track* track2;
    Track* track3;

    // finding kaon (positive) decay vertex
    Int_t vert[nTracks];
    Int_t iLoc = 0;
    for (int it = 0; it < nTracks; it++) {
        track1 = (Track*)branchTrack->At(it);
        if (track1->PID == 321 && Length(track1->X, track1->Y, track1->Z) > 0) {
            vert[iLoc] = it;
            iLoc += 1;
            // cout << "---1" << endl;
        }
    }
    if (iLoc < 1) return iFinalStatesIndexes;

    Int_t foundK = 0;
    Int_t foundPi = 0;
    Int_t iKPos, iKNeg, iPi;
    // finding same vertex kaon and pion
    for (int iar = 0; iar < iLoc; iar++) {
        track1 = (Track*)branchTrack->At(vert[iar]);
        foundK = 0;
        for (int it2 = 0; it2 < nTracks; it2++) {
            // same vertex K (negative)
            track2 = (Track*)branchTrack->At(it2);
            // cout << endl;
            // cout << vert[iar] << "; " << it2 << endl;
            // cout << track1->X << "; " << track1->Y << "; " << track1->Z << "; " << endl;
            // cout << track2->X << "; " << track2->Y << "; " << track2->Z << "; " << endl;
            if (track2->PID == -321 &&
                track1->X == track2->X && track1->Y == track2->Y && track1->Z == track2->Z &&
                track1->Charge + track2->Charge == 0) {
                // cout << "---2" << endl;
                foundPi = 0;
                // same vertex pi
                for (int it3 = 0; it3 < nTracks; it3++) {
                    track3 = (Track*)branchTrack->At(it3);
                    if (abs(track3->PID) == 211 && track1->X == track3->X && track1->Y == track3->Y && track1->Z == track3->Z) {
                        iKPos = vert[iar];
                        iKNeg = it2;
                        iPi = it3;
                        foundK = 1;
                        foundPi = 1;
                        // cout << "---3" << endl;
                        break;
                    }
                }
                if (foundPi == 1) break;
            }
        }
        if (foundPi == 1) break;
    }
    if (not(foundK == 1 && foundPi == 1)) return iFinalStatesIndexes;

    iFinalStatesIndexes.iKPos = iKPos;
    iFinalStatesIndexes.iKNeg = iKNeg;
    iFinalStatesIndexes.iPi = iPi;
    iFinalStatesIndexes.foundFromC = 1;
    return iFinalStatesIndexes;
}

// tag unpaired muon from b-hadron
Int_t findMuIndex(iFinalStates iFS, TClonesArray* branchTrack) {
    Int_t nTracks = branchTrack->GetEntries();
    Track* trackMu;
    Track* trackPi;
    trackPi = (Track*)branchTrack->At(iFS.iPi);
    // find muon
    Int_t foundMu = 0;
    Int_t iMu;
    for (int it = 0; it < nTracks; it++) {
        trackMu = (Track*)branchTrack->At(it);
        if (abs(trackMu->PID) == 13 && trackMu->Charge + trackPi->Charge == 0) {
            iMu = it;
            foundMu += 1;
        }
    }
    if (foundMu != 1) return 99999;
    return iMu;
}

// deduce b-hadron 4 momentum be Z boson 2 body decay
TLorentzVector reconstructB4Momentum(TVector3 v3B, iFinalStates iFS, TClonesArray* branchTrack, TClonesArray* branchEFlowPhoton, TClonesArray* branchEFlowNeutralHadron) {
    const Float_t mK0 = 0.497661;
    const Float_t mBs = 5.36688;
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
                if (it != iFS.iKPos && it != iFS.iKNeg && it != iFS.iPi && it != iFS.iMu) {
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
    Float_t pB = pow(EB * EB - mBs * mBs, 0.5);
    TLorentzVector B;
    Float_t ratio = pB / (pow(v3B.X() * v3B.X() + v3B.Y() * v3B.Y() + v3B.Z() * v3B.Z(), 0.5));
    B.SetPxPyPzE(v3B.X() * ratio, v3B.Y() * ratio, v3B.Z() * ratio, EB);

    return B;
}

// find misID pion, specific to misID bkg type
vector<Int_t> findMisIDPi(iFinalStates iFS, TClonesArray* branchTrack) {
    vector<Int_t> iPiS;
    Int_t nTracks = branchTrack->GetEntries();
    Track* track;
    Track* track0;
    track0 = (Track*)branchTrack->At(iFS.iPi);
    for (int it = 0; it < nTracks; it++) {
        // exclude tagged final states
        if (it == iFS.iKNeg || it == iFS.iKPos || it == iFS.iPi) continue;

        track = (Track*)branchTrack->At(it);
        // find pion
        if (abs(track->PID) != 211) continue;
        // matching charge
        if (track->Charge + track0->Charge != 0) continue;
        // same direction
        if (track->X == 0 && track->Y == 0 && track->Z == 0) continue;

        iPiS.push_back(it);
    }
    return iPiS;
}

