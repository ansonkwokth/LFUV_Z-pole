// ===================
// === R(Lambda_c) ===
// ===================

#include <cmath>
#include <iostream>
#include <random>
#include <vector>
using namespace std;

#include "FeaturesCalculator.C"  // customized feature calculation functions
#include "FeaturesClass.C"       // features list
#include "FinalStatesClass.C"    // final states list
#include "Geometry.C"            // customized geometric functions
#include "Reconstructor.C"       // reconstruction scheme: tagged final states particles and momentum reconstruction related function
#include "Rtypes.h"
#include "SignalBackgroundClassifier.C"  //customized signals and backgrounds classification function, including tagging the truth final state particles
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

    cout << "\n\n\n\n\n\n\n\n\n\n";

    string typeName;
    const char *inputFile;
    const char *outputFile;
    const Float_t noise = noise_ * 0.001;  // changed unit

    // reading data and naming the output file accordingly
    if (type == "s1") {
        cout << "Lambdab->Lambdac tau nu. " << endl;
        inputFile = "./Lambdab_lambdactaunu_pkpi_allchannel_1m_seed0.root";
        if (noise_ == 10) outputFile = "./features/LambdacTauNu_10Noise.root";
        if (noise_ == 20) outputFile = "./features/LambdacTauNu_20Noise.root";

    } else if (type == "s2") {
        cout << "Lambdab->Lambdac mu nu. " << endl;
        inputFile = "./Lambdab_lambdacmunu_pkpi_allchannel_1m_seed0.root";
        if (noise_ == 10) outputFile = "./features/LambdacMuNu_10Noise.root";
        if (noise_ == 20) outputFile = "./features/LambdacMuNu_20Noise.root";

    } else if (type == "b1") {
        cout << "Comb+Cascade Bkg. " << endl;
        inputFile = "./RLambdac_comb_10m_seed1_2_3_4_5_6.root";
        if (noise_ == 10) outputFile = "./features/RLambdacCombCascade_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RLambdacCombCascade_20Noise.root";

    } else if (type == "b2") {
        cout << "Comb Bkg. " << endl;
        inputFile = "./RLambdac_comb_10m_seed1_2_3_4_5_6.root";
        if (noise_ == 10) outputFile = "./features/RLambdacComb_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RLambdacComb_20Noise.root";

    } else if (type == "b3") {
        cout << "Cascade Bkg. " << endl;
        inputFile = "./RLambdac_comb_10m_seed1_2_3_4_5_6.root";
        if (noise_ == 10) outputFile = "./features/RLambdacCascade_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RLambdacCascade_20Noise.root";

    } else if (type == "b4") {
        cout << "Inclusive Bkg. " << endl;
        inputFile = "./RLambdac_comb_10m_seed1_2_3_4_5_6.root";
        if (noise_ == 10) outputFile = "./features/RLambdacInclusive_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RLambdacInclusive_20Noise.root";

    } else if (type == "b5") {
        cout << "MisID Bkg. " << endl;
        inputFile = "./RLambdac_comb_10m_seed1_2_3_4_5_6.root";
        if (noise_ == 10) outputFile = "./features/RLambdacMisID_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RLambdacMisID_20Noise.root";

    } else {
        cout << "Not match. ";
        return;
    }

    // Load lib, and read data
    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    chain.Add(inputFile);
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Int_t numberOfEntries = treeReader->GetEntries();

    // clone the branches
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");                      // truth level
    TClonesArray *branchTrack = treeReader->UseBranch("Track");                            // detector level charged particle
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");                // detector level photon
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");  // detector level neutral hadron

    GenParticle *particle;
    Track *track;
    Tower *eflowphoton;
    Tower *eflowneutralhadron;

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
    Int_t nFoundTracks = 0;
    Int_t nMu = 0;
    Int_t nSameDir = 0;
    Int_t nHcMass = 0;
    Int_t nMuPt = 0;
    Int_t nVert = 0;
    Int_t num = 0;

    // output data
    if (not save) outputFile = "dummy.root";

    TFile fea(outputFile, "recreate");
    TTree tr("t", "features");
    Features *features = new Features;
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
        // progress
        if ((i_en % 100000) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << endl;

        treeReader->ReadEntry(i_en);  // reading the entry

        //==========================================================================
        //===============   Classifying in event type in truth level ===============
        //==========================================================================
        iFinalStates iFSTrue;                  // final states particles, in truth level
        TLorentzVector BTrue, HcTrue, muTrue;  // define lorentz vector for the truth level particle
        TVector3 v3HcTrue, v3MuTrue;
        Int_t passing = 0;
        // check if the event pass the classification of the correspdoning event type
        // and store the truth level info (for signal types)
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

        nEvt += 1;  // count number of truth level events

        //==========================================================================
        //===============   Finding the correspdoning final states   ===============
        //==========================================================================
        // detector level search:
        // finding the correspdoning final states candidates (from c-hadron, and muon)
        // then reconstruct the initial b-hadron

        // finding the indices of those from targeted c-hadron decay first
        iFinalStates iFS = FindFinalStatesIndex(branchTrack);
        if (iFS.foundFromC == 0) continue;  // flag to check all of them are found
        nFoundTracks += 1;

        vector<Int_t> muCandidates;  // indices of muon candidates
        if (type == "b5") {          // if it is misID bkg type, find pions instead of muons
            muCandidates = findMisIDPi(iFS, branchTrack);
        } else {  // for other types, find muons
            muCandidates.push_back(findMuIndex(iFS, branchTrack));
        }

        // loop over muon candidates
        for (Int_t imu : muCandidates) {
            if (imu == 99999) continue;  // if not found any

            iFS.iMu = imu;  // add muon index to the final states list

            //==========================================================================
            //==================   Defining 4 Momentum & 3 Positions  ==================
            //==========================================================================
            Track *pTrack = (Track *)branchTrack->At(iFS.iP);    // proton
            Track *KTrack = (Track *)branchTrack->At(iFS.iK);    // kaon
            Track *piTrack = (Track *)branchTrack->At(iFS.iPi);  // pion
            Track *muTrack = (Track *)branchTrack->At(iFS.iMu);  // muon

            TLorentzVector p;
            p.SetPtEtaPhiM(pTrack->PT, pTrack->Eta, pTrack->Phi, mP);
            TLorentzVector K;
            K.SetPtEtaPhiM(KTrack->PT, KTrack->Eta, KTrack->Phi, mK);
            TLorentzVector pi;
            pi.SetPtEtaPhiM(piTrack->PT, piTrack->Eta, piTrack->Phi, mPi);
            TLorentzVector mu;
            mu.SetPtEtaPhiM(muTrack->PT, muTrack->Eta, muTrack->Phi, mMu);

            // check direction of the tracks
            if (p.Px() * K.Px() + p.Py() * K.Py() + p.Pz() * K.Pz() <= 0) continue;
            if (p.Px() * pi.Px() + p.Py() * pi.Py() + p.Pz() * pi.Pz() <= 0) continue;
            if (p.Px() * mu.Px() + p.Py() * mu.Py() + p.Pz() * mu.Pz() <= 0) continue;
            if (K.Px() * pi.Px() + K.Py() * pi.Py() + K.Pz() * pi.Pz() <= 0) continue;
            if (K.Px() * mu.Px() + K.Py() * mu.Py() + K.Pz() * mu.Pz() <= 0) continue;
            if (pi.Px() * mu.Px() + pi.Py() * mu.Py() + pi.Pz() * mu.Pz() <= 0) continue;
            nSameDir += 1;

            // reconstruct c-hadron
            TLorentzVector Hc;
            Hc = p + K + pi;

            // randomly generate 3D noise vector
            // noise for c-hadron
            Float_t dxC = distribution(genertator);
            Float_t dyC = distribution(genertator);
            Float_t dzC = distribution(genertator);
            TVector3 v3CNoise(dxC, dyC, dzC);
            TVector3 v3C(pTrack->X, pTrack->Y, pTrack->Z);
            v3C += v3CNoise;
            // noise for muon
            Float_t dxMu = distribution(genertator);
            Float_t dyMu = distribution(genertator);
            Float_t dzMu = distribution(genertator);
            TVector3 v3MuNoise(dxMu, dyMu, dzMu);
            TVector3 v3Mu(muTrack->X, muTrack->Y, muTrack->Z);
            v3Mu += v3MuNoise;

            //==========================================================================
            //========================   Vetoing muon and Hc   =========================
            //==========================================================================
            // veto muon
            Float_t disMuTr;
            disMuTr = closestTrack(iFS, mu, v3Mu, noise, branchTrack);
            if (disMuTr < 0.02) continue;
            nMu += 1;

            // veto Hc
            Float_t disHcTr;
            disHcTr = closestTrack(iFS, Hc, v3C, noise, branchTrack);
            if (disHcTr < 0.02) continue;

            //==========================================================================
            //=============================   Apply cuts   =============================
            //==========================================================================
            // c-hadron mass cut
            if (not(Hc.M() >= 2.272 && Hc.M() < 2.3)) continue;
            nHcMass += 1;

            // muon pT cut
            if (not(mu.Pt() > 1.2)) continue;
            nMuPt += 1;

            //==========================================================================
            //=======================   Deduce B decay vertex   ========================
            //==========================================================================
            // c-hadron vertex distance cut
            if (Length(v3C.X(), v3C.Y(), v3C.Z()) < 0.05) continue;

            // tracks distance cut
            Float_t LHcMu, sHc, sMu;
            distance_2lines(v3C.X(), v3C.Y(), v3C.Z(), Hc.Px(), Hc.Py(), Hc.Pz(),
                            v3Mu.X(), v3Mu.Y(), v3Mu.Z(), mu.Px(), mu.Py(), mu.Pz(),
                            &LHcMu, &sHc, &sMu);
            if (LHcMu > 0.5) continue;

            // deduced location (b-hadron decay vertex)
            Float_t X = v3C.X() + sHc * Hc.Px();
            Float_t Y = v3C.Y() + sHc * Hc.Py();
            Float_t Z = v3C.Z() + sHc * Hc.Pz();
            TVector3 v3B(X, Y, Z);
            if (Length(v3C.X(), v3C.Y(), v3C.Z()) < 0.05) continue;
            nVert += 1;

            //==========================================================================
            //===========================   Reconstruct B   ============================
            //==========================================================================
            TLorentzVector B;
            B = reconstructB4Momentum(v3B, iFS, branchTrack, branchEFlowPhoton, branchEFlowNeutralHadron);
            if (isnan(B.P())) continue;
            if (B.Px() == 99999 && B.Py() == 99999 && B.Pz() == 99999 && B.E() == 99999) continue;

            // invariant mass (of reconstructed b-hadron) cut
            TLorentzVector HcMu;
            HcMu = p + K + pi + mu;
            if (HcMu.M() > mLambdab) continue;

            //==========================================================================
            //==============================   Features   ==============================
            //==========================================================================
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
            features->sMinMuBVert = pow(LHcMu, 0.5);
            features->sMinMuHcVert = distance_linepoint(v3C, v3Mu, mu);
            features->sMinMuTr = disMuTr;
            features->sMinHcTr = disHcTr;
            features->sPVHc = Length(v3C.X(), v3C.Y(), v3C.Z());
            features->mHcMu = HcMu.M();

            // p_perp, m_corr
            features->pPerp = cal_pPerp(v3B, Hc, mu);
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

            // truth level info
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

            // K0 info
            TLorentzVector K0S;
            K0S = reconstructK0S(branchTrack, iFS, noise);
            TLorentzVector K0SHcMu;
            K0SHcMu = K0S + Hc + mu;
            features->mK0SHcMu = K0SHcMu.M();
            features->pK0S = K0S.P();

            // distance info
            Float_t sMinMuHcVertTrue;
            sMinMuHcVertTrue = distance_linepoint(v3HcTrue, v3MuTrue, muTrue);
            features->sMinMuHcVertTrue = sMinMuHcVertTrue;

            // fill the output file
            tr.Fill();

            // count the final reconstructed events
            num += 1;
        }
    }

    cout << "\n\n\n\n";
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

