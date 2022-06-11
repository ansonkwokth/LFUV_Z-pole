// ==================
// === R(D_s^(*)) ===
// ==================
//
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

void allinone_storeBVertex(
    const string type,          // data types (signal channel, and background types)
    const Float_t noise_ = 10,  // amount of noise injected to vertex location (unit of microns)
    const Bool_t save = false,  // saving the output features file (.root)
    Int_t num_test = 0,         // number of events to be ran (default 0: run all)
    const Float_t alpha = 1) {  // modified width of \Delta m peak (default 1: not modified)

    cout << "\n\n\n\n\n\n\n\n\n";

    string typeName;
    const char* inputFile;
    const char* outputFile;
    const Float_t noise = noise_ * 0.001;  // change unit

    // reading data and naming the output file accordingly
    if (type == "s1") {
        cout << "Bs->Ds tau nu. " << endl;
        inputFile = "./Bs0Dstaunu-Dsphipi-phiKK_100k_RandomSeed0.root";
        if (noise_ == 10) outputFile = "./features/DsTauNu_10Noise.root";
        if (noise_ == 20) outputFile = "./features/DsTauNu_20Noise.root";

    } else if (type == "s2") {
        cout << "Bs->Ds mu nu. " << endl;
        inputFile = "./Bs0Dsmunu-Dsphipi-phiKK_100k_RandomSeed0.root";
        if (noise_ == 10) outputFile = "./features/DsMuNu_10Noise.root";
        if (noise_ == 20) outputFile = "./features/DsMuNu_20Noise.root";

    } else if (type == "s3") {
        cout << "Bs->Ds* mu nu. " << endl;
        // inputFile = "./Bs0Dsstartaunu-Dsphipi-phiKK_30k_RandomSeed0_30k_RandomSeed1.root";
        inputFile = "./Bs0Dsstartaunu-Dsphipi-phiKK_30k_RandomSeed0_30k_RandomSeed1_30k_RandomSeed2.root";
        if (noise_ == 10) {
            outputFile = "./features/DsstarTauNu_10Noise.root";
            if (abs(alpha - 0) < 1e-6) outputFile = "./features/DsstarTauNu_10Noise_0Alpha.root";
            if (abs(alpha - 0.1) < 1e-6) outputFile = "./features/DsstarTauNu_10Noise_01Alpha.root";
            if (abs(alpha - 0.5) < 1e-6) outputFile = "./features/DsstarTauNu_10Noise_05Alpha.root";
            if (abs(alpha - 2) < 1e-6) outputFile = "./features/DsstarTauNu_10Noise_2Alpha.root";
        } else if (noise_ == 20) {
            outputFile = "./features/DsstarTauNu_20Noise.root";
        }

    } else if (type == "s4") {
        cout << "Bs->Ds* mu nu. " << endl;
        inputFile = "./Bs0Dsstarmunu-Dsphipi-phiKK_30k_RandomSeed0.root";
        if (noise_ == 10) {
            outputFile = "./features/DsstarMuNu_10Noise.root";
            if (abs(alpha - 0) < 1e-6) outputFile = "./features/DsstarMuNu_10Noise_0Alpha.root";
            if (abs(alpha - 0.1) < 1e-6) outputFile = "./features/DsstarMuNu_10Noise_01Alpha.root";
            if (abs(alpha - 0.5) < 1e-6) outputFile = "./features/DsstarMuNu_10Noise_05Alpha.root";
            if (abs(alpha - 2) < 1e-6) outputFile = "./features/DsstarMuNu_10Noise_2Alpha.root";
        } else if (noise_ == 20) {
            outputFile = "./features/DsstarMuNu_20Noise.root";
        }

    } else if (type == "b1") {
        cout << "Comb+Cascade Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        if (noise_ == 10) outputFile = "./features/RDsCombCascade_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RDsCombCascade_20Noise.root";

    } else if (type == "b2") {
        cout << "Comb Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        if (noise_ == 10) outputFile = "./features/RDsComb_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RDsComb_20Noise.root";
        if (noise_ == 0) outputFile = "./features/RDsComb_0Noise.root";

    } else if (type == "b3") {
        cout << "Cascade Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        if (noise_ == 10) outputFile = "./features/RDsCascade_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RDsCascade_20Noise.root";
        if (noise_ == 0) outputFile = "./features/RDsCascade_0Noise.root";

    } else if (type == "b4") {
        cout << "Inclusive Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        if (noise_ == 10) outputFile = "./features/RDsInclusive_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RDsInclusive_20Noise.root";

    } else if (type == "b5") {
        cout << "MisID Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        if (noise_ == 10) outputFile = "./features/RDsMisID_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RDsMisID_20Noise.root";

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
    const Float_t mDsstar = 2.1122;
    const Float_t mBs = 5.36688;

    // Count number of targeted events
    Int_t nEvt = 0;
    Int_t nFoundFinalStates = 0;
    Int_t nRecoHc = 0;
    Int_t nRecoHb = 0;

    Int_t nFoundTracks = 0;
    Int_t nMu = 0;
    Int_t nSameDir = 0;
    Int_t nHcMass = 0;
    Int_t nMuPt = 0;
    Int_t nVert = 0;
    Int_t num = 0;
    Int_t nPi_MisID = 0;  // total number of possible mis-id event (before preselection)

    // output daata
    if (not save) outputFile = "dummy.root";

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
    cout << "# of events:\t" << numberOfEntries << "\n\n";
    if (num_test == 0) num_test = numberOfEntries;

    // loop over events
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // for (Int_t i_en = 1065960; i_en < num_test; i_en++) {
        if ((i_en % 100000) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << endl;

        treeReader->ReadEntry(i_en);  // reading the entry

        //==========================================================================
        //===============   Classifying in event type in truth level ===============
        //==========================================================================
        iFinalStates iFSTrue;                           // final states particles, in truth level
        TLorentzVector BTrue, HcTrue, muTrue, phoTrue;  // define lorentz vectors for the truh level
        TVector3 v3HcTrue, v3MuTrue, v3BTrue;
        Int_t passing = 0;
        Int_t nPi_MisID_i = 0;
        // check if the event pass the classification of the correspdoning event type
        // and store the truth level info (for signal types)
        if (type == "s1") {
            passing = ClassifySignal(branchParticle, &BTrue, &HcTrue, &muTrue, &v3HcTrue, &v3MuTrue, &v3BTrue, &iFSTrue, 1);
        } else if (type == "s2") {
            passing = ClassifySignal(branchParticle, &BTrue, &HcTrue, &muTrue, &v3HcTrue, &v3MuTrue, &v3BTrue, &iFSTrue, 2);
        } else if (type == "s3") {
            passing = ClassifyExcitedSignal(branchParticle, &BTrue, &HcTrue, &muTrue, &phoTrue, &v3HcTrue, &v3MuTrue, &v3BTrue, &iFSTrue, 1);
        } else if (type == "s4") {
            passing = ClassifyExcitedSignal(branchParticle, &BTrue, &HcTrue, &muTrue, &phoTrue, &v3HcTrue, &v3MuTrue, &v3BTrue, &iFSTrue, 2);
        } else if (type == "b1" || type == "b1" || type == "b2" || type == "b3" || type == "b4") {
            passing = ClassifyBkg(branchParticle, type);
        } else if (type == "b5") {
            passing = ClassifyMisID(branchParticle, &nPi_MisID_i);
        } else {
            cout << " No such Signal/Bkg" << endl;
        }
        if (passing == 0) continue;
        // if (not(i_en == 135488 ||
        // i_en == 221766 ||
        // i_en == 418007 ||
        // i_en == 1210 ||
        // i_en == 93994 ||
        // i_en == 814389 ||
        // i_en == 1065964 ||
        // i_en == 152124 ||
        // i_en == 956556 ||
        // i_en == 175734 ||
        // i_en == 850887 ||
        // i_en == 1042933 ||
        // i_en == 15560 ||
        // i_en == 527245 ||
        // i_en == 267106 ||
        // i_en == 1113731 ||
        // i_en == 764191 ||
        // i_en == 991801 ||
        // i_en == 867118 ||
        // i_en == 129808 ||
        // i_en == 716346 ||
        // i_en == 303735 ||
        // i_en == 348217 ||
        // i_en == 728717 ||
        // i_en == 667216 ||
        // i_en == 229967 ||
        // i_en == 876499 ||
        // i_en == 722357 ||
        // i_en == 808232 ||
        // i_en == 357903 ||
        // i_en == 889990 ||
        // i_en == 569199 ||
        // i_en == 1089384 ||
        // i_en == 337216 ||
        // i_en == 1065964 ||
        // i_en == 566768 ||
        // i_en == 95693 ||
        // i_en == 654751 ||
        // i_en == 308169 ||
        // i_en == 1090763 ||
        // i_en == 449866 ||
        // i_en == 110344 ||
        // i_en == 48235 ||
        // i_en == 844066 ||
        // i_en == 180790 ||
        // i_en == 688361 ||
        // i_en == 10221 ||
        // i_en == 770465 ||
        // i_en == 297143 ||
        // i_en == 229946 ||
        // i_en == 503050 ||
        // i_en == 349636 ||
        // i_en == 274912 ||
        // i_en == 190554 ||
        // i_en == 220423 ||
        // i_en == 564127 ||
        // i_en == 38912 ||
        // i_en == 681440 ||
        // i_en == 555113 ||
        // i_en == 577421 ||
        // i_en == 962998 ||
        // i_en == 64512 ||
        // i_en == 880682 ||
        // i_en == 801338 ||
        // i_en == 1064439 ||
        // i_en == 281695 ||
        // i_en == 164195 ||
        // i_en == 165499 ||
        // i_en == 395443 ||
        // i_en == 697040 ||
        // i_en == 470169 ||
        // i_en == 292471)) continue;
        // cout << " .......i_en : " << i_en << "\n";
        //  if (i_en == 1210) cout << "i en: " << i_en << "\n";

        nEvt += 1;  // count number of truth level events
        nPi_MisID += nPi_MisID_i;

        //==========================================================================
        //===============   Finding the correspdoning final states   ===============
        //==========================================================================
        // detector level search:
        // finding the correspdoning final states candidates (from c-hadron, and muon)
        // then reconstruct the initial b-hadron

        // finding the indices of those from targeted c-hadron decay first
        iFinalStates iFS = FindFinalStatesIndex(branchTrack);
        if (iFS.foundFromC == 0) continue;

        vector<Int_t> muCandidates;  // indices of muon candidates
        if (type == "b5") {          // if it is to calculate midIS, find pions instead of muons
            muCandidates = findMisIDPi(iFS, branchTrack);
        } else {  // for other types, find muons
            muCandidates.push_back(findMuIndex(iFS, branchTrack));
        }

        // loop over muon candidates
        for (Int_t imu : muCandidates) {
            if (imu == 99999) continue;  // if not found any

            iFS.iMu = imu;  // add muon index to the final states list
            nFoundFinalStates += 1;

            //==========================================================================
            //==================   Defining 4 Momentum & 3 Positions  ==================
            //==========================================================================
            Track* KNegTrack = (Track*)branchTrack->At(iFS.iKNeg);  // negative kaon
            Track* KPosTrack = (Track*)branchTrack->At(iFS.iKPos);  // positive kaon
            Track* piTrack = (Track*)branchTrack->At(iFS.iPi);      // pion
            Track* muTrack = (Track*)branchTrack->At(iFS.iMu);      // muon

            TLorentzVector KNeg;
            KNeg.SetPtEtaPhiM(KNegTrack->PT, KNegTrack->Eta, KNegTrack->Phi, mK);
            TLorentzVector KPos;
            KPos.SetPtEtaPhiM(KPosTrack->PT, KPosTrack->Eta, KPosTrack->Phi, mK);
            TLorentzVector pi;
            pi.SetPtEtaPhiM(piTrack->PT, piTrack->Eta, piTrack->Phi, mPi);
            TLorentzVector mu;
            mu.SetPtEtaPhiM(muTrack->PT, muTrack->Eta, muTrack->Phi, mMu);

            // reconstruct c-hadron
            TLorentzVector Hc;
            Hc = KNeg + KPos + pi;

            // randomly generate 3D noise vector
            // noise for c-hadron
            Float_t dxC = distribution(genertator);
            Float_t dyC = distribution(genertator);
            Float_t dzC = distribution(genertator);
            TVector3 v3CNoise(dxC, dyC, dzC);
            TVector3 v3C(piTrack->X, piTrack->Y, piTrack->Z);
            v3C += v3CNoise;
            // noise for muon
            Float_t dxMu = distribution(genertator);
            Float_t dyMu = distribution(genertator);
            Float_t dzMu = distribution(genertator);
            TVector3 v3MuNoise(dxMu, dyMu, dzMu);
            TVector3 v3Mu(muTrack->X, muTrack->Y, muTrack->Z);
            if (Length(v3Mu.X(), v3Mu.Y(), v3Mu.Z()) <= 0) continue;  // avoid the PV muon, if adding noise to make it accidently pass the cut
            v3Mu += v3MuNoise;

            // ================================
            // ==========   Ds cut   ==========
            // ================================
            // c-hadron vertex distance cut
            if (Length(v3C.X(), v3C.Y(), v3C.Z()) < 0.5) continue;
            //  KK mass cut
            TLorentzVector KK;
            KK = KPos + KNeg;
            if (not(KK.M() >= 1.008 && KK.M() <= 1.032)) continue;
            //  c-hadron mass cut
            if (not(Hc.M() >= 1.945 && Hc.M() <= 1.995)) continue;
            // veto Hc
            Float_t disHcTr;
            disHcTr = closestTrack(iFS, Hc, v3C, noise, branchTrack);
            if (disHcTr < 0.02) continue;
            nRecoHc += 1;

            // ===========================================
            // ==========   unpaired muon cut   ==========
            // ===========================================
            // need to be from Signal hemisphere
            if (not(v3Mu.X() * v3C.X() + v3Mu.Y() * v3C.Y() + v3Mu.Z() * v3C.Z() > 0)) continue;
            // muon pT cut
            if (not(mu.Pt() > 1.2)) continue;
            // veto muon
            Float_t disMuTr;
            disMuTr = closestTrack(iFS, mu, v3Mu, noise, branchTrack);
            if (disMuTr < 0.02) continue;

            // ===============================================
            // ==========   m(K K pi mu) mass cut   ==========
            // ===============================================
            // invariant mass (of reconstructed b-hadron) cut
            TLorentzVector HcMu;
            HcMu = KNeg + KPos + pi + mu;
            if (HcMu.M() > mBs) continue;

            //==========================================================================
            //=======================   Deduce B decay vertex   ========================
            //==========================================================================
            // tracks distance cut
            Float_t LHcMu, sHc, sMu;
            distance_2lines(v3C.X(), v3C.Y(), v3C.Z(), Hc.Px(), Hc.Py(), Hc.Pz(),
                            v3Mu.X(), v3Mu.Y(), v3Mu.Z(), mu.Px(), mu.Py(), mu.Pz(),
                            &LHcMu, &sHc, &sMu);

            // deduced location
            Float_t X = v3C.X() + sHc * Hc.Px();
            Float_t Y = v3C.Y() + sHc * Hc.Py();
            Float_t Z = v3C.Z() + sHc * Hc.Pz();
            TVector3 v3B(X, Y, Z);

            //==========================================================================
            //===========================   Reconstruct B   ============================
            //==========================================================================
            TLorentzVector B;
            B = reconstructB4Momentum(v3B, iFS, branchTrack, branchEFlowPhoton, branchEFlowNeutralHadron);
            if (isnan(B.P())) continue;
            if (B.Px() == 99999 && B.Py() == 99999 && B.Pz() == 99999 && B.E() == 99999) continue;

            //==========================================================================
            //==============================   Features   ==============================
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
            features->sMinMuBVert = LHcMu;
            features->sMinMuHcVert = distance_linepoint(v3C, v3Mu, mu);
            features->sMinMuTr = disMuTr;
            features->sMinHcTr = disHcTr;
            features->sPVHc = Length(v3C.X(), v3C.Y(), v3C.Z());
            features->mHcMu = HcMu.M();
            // p_perp, m_corr
            // features->pPerpHc = cal_pPerp_Hc(v3B, Hc);
            features->pPerp = cal_pPerp(v3B, Hc, mu);
            features->mCorr = cal_mCorr(v3B, Hc, mu);
            // impact parameters
            impactParams impParams;
            impParams = FindImpactParams(branchTrack, iFS, v3B);
            features->D0Max = impParams.D0Max;
            features->D0Sum = impParams.D0Sum;
            features->DzMax = impParams.DzMax;
            features->DzSum = impParams.DzSum;
            // finding the photon, for modifying the width of \Delta m peak
            whichPhoton wPho;
            wPho = calDeltaM(branchEFlowPhoton, Hc, phoTrue, alpha);
            features->DeltaM = wPho.DeltaM;
            features->correctPhoton = wPho.correctPhoton;
            iFS.iPho = wPho.iPho;
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
            features->EPhoTrue = phoTrue.P();
            Int_t isDsPho;
            if (wPho.iCheatedPho == 99999) {
                isDsPho = 0;
            } else {
                isDsPho = 1;
            }

            features->isDsPho = isDsPho;
            Float_t DeltaRDsPhoTrue = pow(pow(HcTrue.Phi() - phoTrue.Phi(), 2) + pow(HcTrue.Eta() - phoTrue.Eta(), 2), 0.5);
            features->DeltaRDsPhoTrue = DeltaRDsPhoTrue;

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
            Float_t sMinMuBVertTrue = distance_linepoint(v3BTrue, v3MuTrue, muTrue);
            features->sMinMuBVertTrue = sMinMuBVertTrue;
            // cout << " s true: " << sMinMuBVertTrue << "\n";
            // cout << " s : " << features->sMinMuBVert << "\n";
            // cout << endl;

            // fill the output file
            tr.Fill();

            // count the final reconstructed events
            nRecoHb += 1;
        }
    }

    cout << endl;
    cout << endl;
    cout << endl;
    if (type == "b5") {
        cout << "Number of mid-id: \t\t" << nPi_MisID << endl;
    }
    cout << " Number of Events: \t\t" << nEvt << endl;
    cout << " Found all final states: \t" << nFoundFinalStates << endl;
    cout << " Reconstructed Hc: \t\t" << nRecoHc << endl;
    cout << " Number of Reconstructions: \t" << nRecoHb << endl;
    cout << " \t\t\tEff.: \t" << 100 * float(nRecoHb) / float(nEvt) << "%" << endl;

    // cout << "Matched vertex: \t\t" << nFoundTracks << endl;
    // cout << "All in same direction: \t\t" << nSameDir << endl;
    // cout << "Matched muon: \t\t\t" << nMu << endl;
    // cout << "Hc in mass range: \t\t" << nHcMass << endl;
    // cout << "Mu PT cut: \t\t\t" << nMuPt << endl;
    // cout << "Deduced vertex cut: \t\t" << nVert << endl;
    // cout << nFoundFinalStates << endl;

    if (save) {
        tr.Write();
        fea.Close();
    }
    cout << "Writing to:\t" << outputFile << endl;
}

