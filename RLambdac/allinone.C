#include <cmath>
#include <iostream>
#include <random>
#include <vector>
using namespace std;

#include "FeaturesCalculator.h"
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

