
#include <cmath>
#include <iostream>
#include <random>
using namespace std;

#include "Geometry.h"
#include "Rtypes.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

Float_t cal_pPerp(TVector3 v3B, TLorentzVector Hc, TLorentzVector Mu) {
    TLorentzVector HcMu;
    HcMu = Hc + Mu;
    Float_t proj = (HcMu.Px() * v3B.X() + HcMu.Py() * v3B.Y() + HcMu.Pz() * v3B.Z()) / pow(Length(v3B.X(), v3B.Y(), v3B.Z()), 2);
    Float_t pPerp = Length(HcMu.Px() - proj * v3B.X(), HcMu.Py() - proj * v3B.Y(), HcMu.Pz() - proj * v3B.Z());
    return pPerp;
}

Float_t cal_mCorr(TVector3 v3B, TLorentzVector Hc, TLorentzVector Mu) {
    Float_t pPerp = cal_pPerp(v3B, Hc, Mu);
    TLorentzVector HcMu;
    HcMu = Hc + Mu;
    Float_t mCorr = pow(HcMu.M() * HcMu.M() + pPerp * pPerp, 0.5) + pPerp;
    return mCorr;
}

struct impactParams {
    Float_t D0Max = 0;
    Float_t D0Sum = 0;
    Float_t DzMax = 0;
    Float_t DzSum = 0;
};

impactParams FindImpactParams(TClonesArray *branchTrack, iFinalStates iFS, TVector3 v3B) {
    impactParams impParams;

    Int_t nTracks = branchTrack->GetEntries();
    Track *track;
    for (int it = 0; it < nTracks; it++) {
        track = (Track *)branchTrack->At(it);
        if (it == iFS.iP || it == iFS.iK || it == iFS.iPi || it == iFS.iMu) {
            continue;
        }
        TLorentzVector impTr;
        impTr.SetPtEtaPhiE(track->PT, track->Eta, track->Phi, track->P);
        if ((impTr.Px() * v3B.X() + impTr.Py() * v3B.Y() + impTr.Pz() * v3B.Z() > 0) &&
            (track->X * v3B.X() + track->Y * v3B.Y() + track->Z * v3B.Z() > 0)) {
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

isolationVars FindIsolationVars(
    TClonesArray *branchTrack,
    TClonesArray *branchEFlowNeutralHadron,
    TClonesArray *branchEFlowPhoton,
    iFinalStates iFS,
    TLorentzVector Hc) {
    isolationVars isoVars;

    const Float_t mK0 = 0.497661;
    TLorentzVector B03Neu, B06Neu, B03NeuHad, B06NeuHad, B03NeuPho, B06NeuPho;
    TLorentzVector B03Chg, B06Chg, B03ChgPV, B06ChgPV, B03ChgDV, B06ChgDV;

    Int_t nEFlowNeutralHadrons = branchEFlowNeutralHadron->GetEntries();
    Int_t nEFlowPhotons = branchEFlowPhoton->GetEntries();
    Int_t nTracks = branchTrack->GetEntries();

    Tower *eflowNeuHad;
    Tower *eflowPho;
    Track *track;

    for (int in = 0; in < nEFlowNeutralHadrons; in++) {
        eflowNeuHad = (Tower *)branchEFlowNeutralHadron->At(in);
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
        eflowPho = (Tower *)branchEFlowPhoton->At(iph);
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
        if (it == iFS.iP || it == iFS.iK || it == iFS.iPi || it == iFS.iMu) {
            continue;
        }
        track = (Track *)branchTrack->At(it);
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

TLorentzVector reconstructK0S(TClonesArray *branchTrack, iFinalStates iFS, double noise) {
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
    Track *track;
    Track *trackPos;
    Track *trackNeg;
    for (int it = 0; it < nTracks; it++) {
        if (it == iFS.iP || it == iFS.iK || it == iFS.iP) {
            continue;
        }
        track = (Track *)branchTrack->At(it);
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
        trackPos = (Track *)branchTrack->At(ipos);
        for (int ineg : iNegS) {
            trackNeg = (Track *)branchTrack->At(ineg);
            if (not(trackPos->X == trackNeg->X && trackPos->Y && trackNeg->Y && trackPos->Z == trackNeg->Z)) {
                continue;
            }

            Float_t XK0S = trackPos->X + distribution(genertator);
            Float_t YK0S = trackPos->Y + distribution(genertator);
            Float_t ZK0S = trackPos->Z + distribution(genertator);

            if (Length(XK0S, YK0S, ZK0S) <= 10) {
                continue;
            }
            // cout << "Pos pi vertex length: " << Length(trackPos->X, trackPos->Y, trackPos->Z) << endl;
            // cout << "Neg pi vertex length: " << Length(trackNeg->X, trackNeg->Y, trackNeg->Z) << endl;
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
    Track *trackOther;
    TLorentzVector target = K0S;
    TVector3 v3Target = v3K0S;

    for (int it = 0; it < nTracks; it++) {
        if (it == iFS.iP || it == iFS.iK || it == iFS.iPi || it == iFS.iMu || it == savediNeg || it == savediPos) {
            continue;
        }
        trackOther = (Track *)branchTrack->At(it);
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
        distance_2lines(XTr, YTr, ZTr, otherTrack.Px(), otherTrack.Py(), otherTrack.Pz(),
                        v3Target.X(), v3Target.Y(), v3Target.Z(), target.Px(), target.Py(), target.Pz(),
                        &L, &s1, &s2);
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

