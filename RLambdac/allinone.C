#include <cmath>
#include <iostream>
#include <random>
#include <vector>
using namespace std;

#include "FeaturesClass.C"
#include "FinalStatesClass.C"
#include "Geometry.h"
#include "Reconstructor.h"
#include "Rtypes.h"
#include "SignalBackgroundClassifier.h"
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
#include <vector>

// Load Delphes for reading the data
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

//---------- {{{
Float_t cal_pPerp(TVector3 v3B, TLorentzVector Hc, TLorentzVector Mu) {
    TLorentzVector HcMu;
    HcMu = Hc + Mu;
    Float_t proj = (HcMu.Px() * v3B.X() + HcMu.Py() * v3B.Y() + HcMu.Pz() * v3B.Z()) / pow(Length(v3B.X(), v3B.Y(), v3B.Z()), 2);
    Float_t pPerp = Length(HcMu.Px() - proj * v3B.X(), HcMu.Py() - proj * v3B.Y(), HcMu.Pz() - proj * v3B.Z());
    return pPerp;
}
//}}}

//---------- {{{
Float_t cal_mCorr(TVector3 v3B, TLorentzVector Hc, TLorentzVector Mu) {
    Float_t pPerp = cal_pPerp(v3B, Hc, Mu);
    TLorentzVector HcMu;
    HcMu = Hc + Mu;
    Float_t mCorr = pow(HcMu.M() * HcMu.M() + pPerp * pPerp, 0.5) + pPerp;
    return mCorr;
}
//}}}

//---------- {{{
struct impactParams {
    Float_t D0Max = 0;
    Float_t D0Sum = 0;
    Float_t DzMax = 0;
    Float_t DzSum = 0;
};
//}}}

//---------- {{{
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
//}}}

//---------- {{{
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
Float_t Angle(
    Float_t X1, Float_t Y1, Float_t Z1,
    Float_t X2, Float_t Y2, Float_t Z2) {
    Float_t angle = acos((X1 * X2 + Y1 * Y2 + Z1 * Z2) / (Length(X1, Y1, Z1) * Length(X2, Y2, Z2)));
    return angle;
}
//}}}

//---------- {{{
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
//}}}

// {{{
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
//}}}

void allinone(
    const string type,
    const Float_t noise_ = 10,
    const Bool_t save = true,
    Int_t num_test = 0) {
    cout << "\n\n\n\n\n\n\n\n\n\n"
         << endl;
    string typeName;
    const char *inputFile;
    const char *outputFile;
    const Float_t noise = noise_ * 0.001;
    if (type == "s1") {
        typeName = "Lambdab->Lambdac tau nu. ";
        cout << "Lambdab->Lambdac tau nu. " << endl;
        inputFile = "./Lambdab_lambdactaunu_pkpi_allchannel_1m_seed0.root";
        if (noise_ == 10) {
            outputFile = "./features/LambdacTauNu_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/LambdacTauNu_20Noise.root";
        }
    } else if (type == "s2") {
        cout << "Lambdab->Lambdac mu nu. " << endl;
        inputFile = "./Lambdab_lambdacmunu_pkpi_allchannel_1m_seed0.root";
        if (noise_ == 10) {
            outputFile = "./features/LambdacMuNu_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/LambdacMuNu_20Noise.root";
        }
    } else if (type == "b1") {
        cout << "Comb+Cascade Bkg. " << endl;
        inputFile = "./RLambdac_comb_10m_seed1_2_3_4_5_6.root";
        outputFile = "./testoutput.root";
        if (noise_ == 10) {
            outputFile = "./features/RLambdacCombCascade_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/RLambdacCombCascade_20Noise.root";
        }
    } else if (type == "b2") {
        cout << "Comb Bkg. " << endl;
        inputFile = "./RLambdac_comb_10m_seed1_2_3_4_5_6.root";
        outputFile = "./testoutput.root";
        if (noise_ == 10) {
            outputFile = "./features/RLambdacComb_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/RLambdacComb_20Noise.root";
        }
    } else if (type == "b3") {
        cout << "Cascade Bkg. " << endl;
        inputFile = "./RLambdac_comb_10m_seed1_2_3_4_5_6.root";
        outputFile = "./testoutput.root";
        if (noise_ == 10) {
            outputFile = "./features/RLambdacCascade_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/RLambdacCascade_20Noise.root";
        }
    } else if (type == "b4") {
        cout << "Inclusive Bkg. " << endl;
        inputFile = "./RLambdac_comb_10m_seed1_2_3_4_5_6.root";
        outputFile = "./testoutput.root";
        if (noise_ == 10) {
            outputFile = "./features/RLambdacInclusive_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/RLambdacInclusive_20Noise.root";
        }
    } else if (type == "b5") {
        cout << "MisID Bkg. " << endl;
        inputFile = "./RLambdac_comb_10m_seed1_2_3_4_5_6.root";
        outputFile = "./testoutput.root";
        if (noise_ == 10) {
            outputFile = "./features/RLambdacMisID_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/RLambdacMisID_20Noise.root";
        }
    } else {
        cout << "Not match. ";
    }

    // Load lib, and read data
    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    chain.Add(inputFile);
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Int_t numberOfEntries = treeReader->GetEntries();

    // clone the branches
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

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

    GenParticle *particle;
    Track *track;
    Tower *eflowphoton;
    Tower *eflowneutralhadron;

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
    Features *features = new Features;
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

    // for (Int_t i_en=0; i_en< numberOfEntries; i_en++) {
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // if (i_en % 1000 == 0) {cout << "\tReconstruction Progress: " << i_en << "/" << numberOfEntries; }
        // if ((i_en % 1000) == 0) {cout << "\rReconstruction Progress: " << i_en << "/" << numberOfEntries; }
        if ((i_en % 100000) == 0) {
            cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << endl;
        }
        // if (i_en >= num_test) { break; }
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

        if (passing == 0) {
            continue;
        }

        // cout << "\n===========================================Event: " << i_en << endl;
        nEvt += 1;
        //==========================================================================
        //===============   Finding the correspdoning final states   ===============
        //==========================================================================
        iFinalStates iFS = FindFinalStatesIndex(branchTrack);
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
        for (Int_t imu : muCandidates) {
            if (imu == 99999) {
                continue;
            }
            // cout << iloop << endl;
            iFS.iMu = imu;
            // cout << " ..." << iFS.iP << "; " << iFS.iK << "; " << iFS.iPi << "; " << iFS.iMu << endl;

            //==========================================================================
            //==================   Defining 4 Momentum & 3 Positions  ==================
            //==========================================================================
            Track *pTrack = (Track *)branchTrack->At(iFS.iP);
            Track *KTrack = (Track *)branchTrack->At(iFS.iK);
            Track *piTrack = (Track *)branchTrack->At(iFS.iPi);
            Track *muTrack = (Track *)branchTrack->At(iFS.iMu);

            TLorentzVector p;
            p.SetPtEtaPhiM(pTrack->PT, pTrack->Eta, pTrack->Phi, mP);
            TLorentzVector K;
            K.SetPtEtaPhiM(KTrack->PT, KTrack->Eta, KTrack->Phi, mK);
            TLorentzVector pi;
            pi.SetPtEtaPhiM(piTrack->PT, piTrack->Eta, piTrack->Phi, mPi);
            TLorentzVector mu;
            mu.SetPtEtaPhiM(muTrack->PT, muTrack->Eta, muTrack->Phi, mMu);

            if (p.Px() * K.Px() + p.Py() * K.Py() + p.Pz() * K.Pz() <= 0) {
                continue;
            }
            if (p.Px() * pi.Px() + p.Py() * pi.Py() + p.Pz() * pi.Pz() <= 0) {
                continue;
            }
            if (p.Px() * mu.Px() + p.Py() * mu.Py() + p.Pz() * mu.Pz() <= 0) {
                continue;
            }
            if (K.Px() * pi.Px() + K.Py() * pi.Py() + K.Pz() * pi.Pz() <= 0) {
                continue;
            }
            if (K.Px() * mu.Px() + K.Py() * mu.Py() + K.Pz() * mu.Pz() <= 0) {
                continue;
            }
            if (pi.Px() * mu.Px() + pi.Py() * mu.Py() + pi.Pz() * mu.Pz() <= 0) {
                continue;
            }

            nSameDir += 1;

            TLorentzVector Hc;
            Hc = p + K + pi;

            Float_t dxC = distribution(genertator);
            Float_t dyC = distribution(genertator);
            Float_t dzC = distribution(genertator);
            Float_t dxMu = distribution(genertator);
            Float_t dyMu = distribution(genertator);
            Float_t dzMu = distribution(genertator);
            TVector3 v3CNoise(dxC, dyC, dzC);
            TVector3 v3MuNoise(dxMu, dyMu, dzMu);
            TVector3 v3C(pTrack->X, pTrack->Y, pTrack->Z);
            TVector3 v3Mu(muTrack->X, muTrack->Y, muTrack->Z);
            v3Mu += v3MuNoise;
            v3C += v3CNoise;

            //==========================================================================
            //========================   Vetoing muon and Hc   =========================
            //==========================================================================

            // veto muon
            Float_t disMuTr;
            disMuTr = closestTrack(iFS, mu, v3Mu, noise, branchTrack);
            if (disMuTr < 0.02) {
                continue;
            }
            nMu += 1;

            // veto Hc
            Float_t disHcTr;
            disHcTr = closestTrack(iFS, Hc, v3C, noise, branchTrack);
            if (disHcTr < 0.02) {
                continue;
            }

            //==========================================================================
            //=============================   Apply cuts   =============================
            //==========================================================================
            if (not(Hc.M() >= 2.272 && Hc.M() < 2.3)) {
                continue;
            }
            nHcMass += 1;
            if (not(mu.Pt() > 1.2)) {
                continue;
            }
            nMuPt += 1;

            //==========================================================================
            //=======================   Deduce B decay vertex   ========================
            //==========================================================================
            if (Length(v3C.X(), v3C.Y(), v3C.Z()) < 0.05) {
                continue;
            }
            Float_t LHcMu, sHc, sMu;

            distance_2lines(v3C.X(), v3C.Y(), v3C.Z(), Hc.Px(), Hc.Py(), Hc.Pz(),
                            v3Mu.X(), v3Mu.Y(), v3Mu.Z(), mu.Px(), mu.Py(), mu.Pz(),
                            &LHcMu, &sHc, &sMu);
            if (LHcMu > 0.5) {
                continue;
            }
            // if (LHcMu > 0.1) { continue; }

            // deduced location
            Float_t X = v3C.X() + sHc * Hc.Px();
            Float_t Y = v3C.Y() + sHc * Hc.Py();
            Float_t Z = v3C.Z() + sHc * Hc.Pz();
            TVector3 v3B(X, Y, Z);

            if (Length(v3C.X(), v3C.Y(), v3C.Z()) < 0.05) {
                continue;
            }
            nVert += 1;

            TLorentzVector B;
            B = reconstructHc4Momentum(v3B, iFS, branchTrack, branchEFlowPhoton, branchEFlowNeutralHadron);
            if (isnan(B.P())) {
                continue;
            }
            if (B.Px() == 99999 && B.Py() == 99999 && B.Pz() == 99999 && B.E() == 99999) {
                continue;
            }

            TLorentzVector HcMu;
            HcMu = p + K + pi + mu;
            if (HcMu.M() > mLambdab) {
                continue;
            }

            //==========================================================================
            //==============================   Features   ==============================
            //==========================================================================
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
            features->sMinMuBVert = pow(LHcMu, 0.5);
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
            if (q2True < 0) {
                cout << "" << endl;
                cout << i_en << "; " << q2True << "; " << features->q2 << endl;
                cout << muTrue.Px() << "; " << muTrue.Py() << "; " << muTrue.Pz() << endl;
                cout << mu.Px() << "; " << mu.Py() << "; " << mu.Pz() << endl;

                cout << HcTrue.Px() << "; " << HcTrue.Py() << "; " << HcTrue.Pz() << endl;
                cout << Hc.Px() << "; " << Hc.Py() << "; " << Hc.Pz() << endl;

                cout << BTrue.Px() << "; " << BTrue.Py() << "; " << BTrue.Pz() << endl;
                cout << B.Px() << "; " << B.Py() << "; " << B.Pz() << endl;
            }
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
            num += 1;
        }
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
    }
    cout << "Writing to:\t" << outputFile << endl;
}

