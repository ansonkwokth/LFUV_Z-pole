#pragma once
#include <cmath>
#include <random>

#include "FinalStatesClass.C"
#include "Rtypes.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

//---------- calculating the distance from the origin
Float_t Length(Float_t X, Float_t Y, Float_t Z) {
    return pow(X * X + Y * Y + Z * Z, 0.5);
}

void distance_2lines(Float_t X1, Float_t Y1, Float_t Z1, Float_t Px1, Float_t Py1, Float_t Pz1, Float_t X2, Float_t Y2, Float_t Z2, Float_t Px2, Float_t Py2, Float_t Pz2,
                     Float_t *L, Float_t *s1, Float_t *s2) {
    /*calculate the closest distance between 2 lines*/

    Float_t dx = X1 - X2;
    Float_t dy = Y1 - Y2;
    Float_t dz = Z1 - Z2;
    // momemtum magnitude
    Float_t P1 = pow(Px1 * Px1 + Py1 * Py1 + Pz1 * Pz1, 0.5);
    Float_t P2 = pow(Px2 * Px2 + Py2 * Py2 + Pz2 * Pz2, 0.5);

    Float_t A = (dx * Px2 + dy * Py2 + dz * Pz2) / (P2 * P2) * (Px1 * Px2 + Py1 * Py2 + Pz1 * Pz2);
    Float_t B = (dx * Px1 + dy * Py1 + dz * Pz1);
    Float_t C = pow((Px1 * Px2 + Py1 * Py2 + Pz1 * Pz2), 2) / (P1 * P1 * P2 * P2);

    // parametrized 2 tracks
    Float_t s1_ = 1 / (P1 * P1) * (A - B) / (1 - C);
    Float_t s2_ = (B + P1 * P1 * s1_) / (Px1 * Px2 + Py1 * Py2 + Pz1 * Pz2);

    // distance squared between 2 tracks
    Float_t L_ = pow((dx + s1_ * Px1 - s2_ * Px2), 2) + pow((dy + s1_ * Py1 - s2_ * Py2), 2) + pow((dz + s1_ * Pz1 - s2_ * Pz2), 2);
    // cout << L_ << "; " << s1_ << "; " << s2_ << endl;

    *L = pow(L_, 0.5);  // change length squred to length
    *s1 = s1_;
    *s2 = s2_;
}

Float_t distance_linepoint(TVector3 v3Point, TVector3 v3Line, TLorentzVector line) {
    Float_t proj = (line.Px() * (v3Point.X() - v3Line.X()) + line.Py() * (v3Point.Y() - v3Line.Y()) + line.Pz() * (v3Point.Z() - v3Line.Z())) / (line.P() * line.P());
    Float_t X = proj * line.Px() - (v3Point.X() - v3Line.X());
    Float_t Y = proj * line.Py() - (v3Point.Y() - v3Line.Y());
    Float_t Z = proj * line.Pz() - (v3Point.Z() - v3Line.Z());
    Float_t s = Length(X, Y, Z);
    return s;
}

//---------- used in veto (used in muon and Hc) {{{
Float_t closestTrack(
    iFinalStates iFS,
    TLorentzVector target,
    TVector3 v3Target,
    double noise,
    TClonesArray *branchTrack) {
    Float_t disTargetTr = 99999;  // closest distance between a track and the target
    Track *trackOther;

    // inject noise, in track location
    std::default_random_engine genertator;
    std::normal_distribution<double> distribution(0, noise / pow(3, 0.5));

    Int_t nTracks = branchTrack->GetEntries();
    for (int it = 0; it < nTracks; it++) {
        if (it == iFS.iP || it == iFS.iK || it == iFS.iPi || it == iFS.iMu) {
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
    return disTargetTr;
}

Float_t Angle(
    Float_t X1, Float_t Y1, Float_t Z1,
    Float_t X2, Float_t Y2, Float_t Z2) {
    Float_t angle = acos((X1 * X2 + Y1 * Y2 + Z1 * Z2) / (Length(X1, Y1, Z1) * Length(X2, Y2, Z2)));
    return angle;
}
