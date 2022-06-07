// =================
// === R(J/\psi) ===
// =================

#include <cmath>
#include <iostream>
#include <random>
#include <vector>
using namespace std;

#include "FeaturesCalculator.C"
#include "FeaturesClass.C"
#include "FinalStatesClass.C"
#include "Geometry.C"
#include "Reconstructor.C"
#include "Rtypes.h"
#include "SignalBackgroundClassifier.C"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

R__ADD_LIBRARY_PATH($DELPHES)
R__LOAD_LIBRARY(libDelphes)

void allinone(
    const string type,          // data types (signal channel, and background types)
    const Float_t noise_ = 10,  // amount of noise injected to vertex location (unit of microns)
    const Bool_t save = true,   // saving the output features file (.root)
    Int_t num_test = 0) {       // number of events to be ran (default 0: run all)

    cout << "\n\n\n\n\n\n\n\n\n\n\n";

    string typeName;
    const char* inputFile;
    const char* outputFile;
    const Float_t noise = noise_ * 0.001;  // change unit

    if (type == "s1") {
        cout << "Bc->Jpsi tau nu. " << endl;
        inputFile = "./Bcjpsitaunu_50m.root";
        if (noise_ == 10) outputFile = "./features/JpsiTauNu_10Noise_NoVeto.root";
        if (noise_ == 20) outputFile = "./features/JpsiTauNu_20Noise_NoVeto.root";

    } else if (type == "s2") {
        cout << "Bc->Jpsi mu nu. " << endl;
        inputFile = "./BcJpsimunu0-2.root";
        if (noise_ == 10) outputFile = "./features/JpsiMuNu_10Noise_NoVeto.root";
        if (noise_ == 20) outputFile = "./features/JpsiMuNu_20Noise_NoVeto.root";

    } else if (type == "b1") {
        cout << "Comb+Cascade Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        if (noise_ == 10) outputFile = "./features/RJpsiCombCascade_10Noise_NoVeto.root";
        if (noise_ == 20) outputFile = "./features/RJpsiCombCascade_20Noise_NoVeto.root";

    } else if (type == "b2") {
        cout << "Comb Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        if (noise_ == 10) outputFile = "./features/RJpsiComb_10Noise_NoVeto.root";
        if (noise_ == 20) outputFile = "./features/RJpsiComb_20Noise_NoVeto.root";
        if (noise_ == 0) outputFile = "./features/RJpsiComb_0Noise_NoVeto.root";

    } else if (type == "b3") {
        cout << "Cascade Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        if (noise_ == 10) outputFile = "./features/RJpsiCascade_10Noise_NoVeto.root";
        if (noise_ == 20) outputFile = "./features/RJpsiCascade_20Noise_NoVeto.root";
        if (noise_ == 0) outputFile = "./features/RJpsiCascade_0Noise_NoVeto.root";

    } else if (type == "b4") {
        cout << "Inclusive Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        // inputFile = "./RJpsi_comb_200m_seed2.root";
        if (noise_ == 10) outputFile = "./features/RJpsiInclusive_10Noise_NoVeto_seed1.root";
        // if (noise_ == 10) outputFile = "./features/RJpsiInclusive_10Noise_NoVeto_seed2.root";
        if (noise_ == 20) outputFile = "./features/RJpsiInclusive_20Noise_NoVeto_seed1.root";
        // if (noise_ == 20) outputFile = "./features/RJpsiInclusive_20Noise_NoVeto_seed2.root";

    } else if (type == "b5") {
        cout << "MisID Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        if (noise_ == 10) outputFile = "./features/RJpsiMisID_10Noise_NoVeto.root";
        if (noise_ == 20) outputFile = "./features/RJpsiMisID_20Noise_NoVeto.root";

    } else {
        cout << "Not match. ";
        return;
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

    GenParticle* particle;
    Track* track;
    Tower* eflowphoton;
    Tower* eflowneutralhadron;

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

    // Count number of targeted events
    Int_t nEvt = 0;
    Int_t nFoundFinalStates = 0;
    Int_t nRecoHc = 0;
    Int_t nRecoHb = 0;
    Int_t nPi_MisID = 0;  // total number of possible mis-id event (before preselection)

    // output daata
    if (not save) outputFile = "dummy.root";

    TFile fea(outputFile, "recreate");
    TTree tr("t", "features");
    Features* features = new Features;
    tr.Branch("features", &features);

    cout << "\nReconstructing: " << typeName << endl;
    cout << "Reading from:\t" << inputFile << endl;
    cout << "Injected noise:\t" << noise_ << " microns (" << noise << "mm)" << endl;
    cout << "Save: " << save << endl;
    cout << "Writing to:\t" << outputFile << endl;
    cout << "# of events:\t" << numberOfEntries << "\n\n";
    if (num_test == 0) num_test = numberOfEntries;

    // loop over events
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // for (Int_t i_en = 4158740; i_en < num_test; i_en++) {
        if ((i_en % 100000) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << endl;

        // if (i_en != 4162471 && i_en != 4171923 && i_en != 4181737 && i_en != 4182789 && i_en != 4184360 && i_en != 4214889 && i_en != 4218443 && i_en != 4225089) continue;
        // if (i_en != 4159035 && i_en != 4159495) continue;
        // cout << " ....: ";
        treeReader->ReadEntry(i_en);  // reading the entry

        //==========================================================================
        //===============   Classifying in event type in truth level ===============
        //==========================================================================
        iFinalStates iFSTrue;                  // final states, in truth level
        TLorentzVector BTrue, HcTrue, muTrue;  // define lorentz vector for the truth level
        TVector3 v3HcTrue, v3MuTrue;
        Int_t passing = 0;
        Int_t nPi_MisID_i = 0;

        if (type == "s1") {
            passing = ClassifySignal(branchParticle, &BTrue, &HcTrue, &muTrue, &v3HcTrue, &v3MuTrue, &iFSTrue, 1);
        } else if (type == "s2") {
            passing = ClassifySignal(branchParticle, &BTrue, &HcTrue, &muTrue, &v3HcTrue, &v3MuTrue, &iFSTrue, 2);
        } else if (type == "b1" || type == "b1" || type == "b2" || type == "b3" || type == "b4") {
            passing = ClassifyBkg(branchParticle, type);
        } else if (type == "b5") {
            passing = ClassifyMisID(branchParticle, &nPi_MisID_i);
        } else {
            cout << " No such Signal/Bkg" << endl;
        }

        if (passing == 0) continue;
        // if (not(i_en == 160734 ||
        // i_en == 2579758 ||
        // i_en == 2192763 ||
        // i_en == 761722 ||
        // i_en == 2447002 ||
        // i_en == 1410482 ||
        // i_en == 3517198 ||
        // i_en == 4158741 ||
        // i_en == 2690601)) continue;
        // cout << " i_en ..............: " << i_en << "\n";

        nEvt += 1;
        nPi_MisID += nPi_MisID_i;

        //==========================================================================
        //===============   Finding the correspdoning final states   ===============
        //==========================================================================
        Int_t nMu = 3;
        if (type == "b5") nMu = 2;  // misID background
        iFinalStates iFS = FindFinalStatesIndex(branchTrack, nMu = nMu);
        if (iFS.foundFromC == 0) continue;

        vector<Int_t> muCandidates;
        if (type == "b5") {
            muCandidates = findMisIDPi(iFS, branchTrack);
        } else {
            muCandidates.push_back(findMuIndex(iFS, branchTrack));
        }

        for (Int_t imu : muCandidates) {
            // no unpaired muon found in the previous step
            if (imu == 99999) continue;

            iFS.iMu = imu;
            nFoundFinalStates += 1;

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
            if (Length(v3Mu.X(), v3Mu.Y(), v3Mu.Z()) <= 0) continue;  // avoid the PV muon, if adding noise to make it accidently pass the cut
            v3Mu += v3MuNoise;
            v3C += v3CNoise;

            // ===================================
            // ==========   J/psi cut   ==========
            // ===================================
            // distance between PV cut
            if (Length(v3C.X(), v3C.Y(), v3C.Z()) < 0.1) continue;
            // momentum cut of the muon, Hc, and mass cut of Hc
            if (not(abs(Hc.M() - 3.09) < 0.055 / 2 && Hc.Pt() > 2.0 / 2 && muPos.P() > 5.0 / 2 && muNeg.P() > 5.0 / 2 && (muPos.Pt() > 1.5 / 2 || muNeg.Pt() > 1.5 / 2))) continue;
            nRecoHc += 1;

            // ===========================================
            // ==========   unpaired muon cut   ==========
            // ===========================================
            // testing
            // if (mu.E() < 2) cout << "Muon E: " << mu.E() << "\n";
            // if (mu.Pt() < 0.1) cout << "Muon P: " << mu.Pt() << "\n";
            //  end testing
            //  need to be from Signal hemisphere
            if (not(v3Mu.X() * v3C.X() + v3Mu.Y() * v3C.Y() + v3Mu.Z() * v3C.Z() > 0)) continue;
            // momentum cut
            if (not(mu.Pt() > 0.75 / 2 && mu.P() > 3.0 / 2)) continue;

            // =========================================
            // ==========   m(3mu) mass cut   ==========
            // =========================================
            TLorentzVector HcMu;
            HcMu = muPos + muNeg + mu;
            if (HcMu.M() > mBc) continue;

            //=====================================================================
            //=======================   Reconstructing B   ========================
            //=====================================================================
            // deduced location (In J/psi case, B decay vertex = Jpsi decay vertex)
            Float_t X = v3C.X();
            Float_t Y = v3C.Y();
            Float_t Z = v3C.Z();
            TVector3 v3B(X, Y, Z);

            TLorentzVector B;
            B = reconstructB4Momentum(v3B, iFS, branchTrack, branchEFlowPhoton, branchEFlowNeutralHadron);
            if (isnan(B.P())) continue;
            if (B.Px() == 99999 && B.Py() == 99999 && B.Pz() == 99999 && B.E() == 99999) continue;

            //==========================================================================
            //==============================   Feature    ==============================
            //==========================================================================
            // event index
            features->iEvt = i_en;
            // q2
            TLorentzVector q;
            q = B - Hc;
            Float_t q2 = q.E() * q.E() - q.P() * q.P();
            if (q2 < -10 || q2 > 15) continue;
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

            // veto muon
            Float_t disMuTr;
            disMuTr = closestTrack(iFS, mu, v3Mu, noise, branchTrack);
            // veto Hc
            Float_t disHcTr;
            disHcTr = closestTrack(iFS, Hc, v3C, noise, branchTrack);
            features->sMinMuTr = disMuTr;
            features->sMinHcTr = disHcTr;
            features->sPVHc = Length(v3C.X(), v3C.Y(), v3C.Z());
            features->mHcMu = HcMu.M();
            // p_perp
            // features->pPerpHc = cal_pPerp_Hc(v3B, Hc);
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
            nRecoHb += 1;
        }
    }

    cout << endl;
    cout << endl;
    cout << endl;
    cout << "\nReconstructing: " << typeName << endl;
    if (type == "b5") {
        cout << " Number of mid-id: \t\t" << nPi_MisID << endl;
    }
    cout << " Number of Events: \t\t" << nEvt << endl;
    cout << " Found all final states: \t" << nFoundFinalStates << endl;
    cout << " Reconstructed Hc: \t\t" << nRecoHc << endl;
    cout << " Number of Reconstructions: \t" << nRecoHb << endl;
    cout << " \t\t\tEff.: \t" << 100 * float(nRecoHb) / float(nEvt) << "%" << endl;
    // cout << " Mu PT cut: \t\t\t" << nMuPt << endl;
    // cout << " Deduced vertex cut: \t\t" << nVert << endl;
    // cout << " Number of Reconstrcutions: \t" << num << " / " << nEvt << "\t" << nHcMass << " > " << num << "; " << float(num) / float(nHcMass) * 100 << endl;

    if (save) {
        tr.Write();
        fea.Close();
        cout << "Writing to:\t" << outputFile << endl;
    }
}

