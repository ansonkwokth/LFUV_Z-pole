#include <cmath>
#include <iostream>
using namespace std;

#include "FinalStatesClass.C"
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

//{{{
Float_t calCosinceTheta(TLorentzVector p1, TLorentzVector p2) {
    return (p1.Px() * p2.Px() + p1.Py() * p2.Py() + p1.Pz() * p2.Pz()) / (p1.P() * p2.P());
}
//}}}

// calculate number of b's and c's{{{
void calBsCs(TClonesArray* branchParticle,
             Int_t* nB_, Int_t* nC_,
             Float_t* E1_, Float_t* E2_,
             Float_t* angle1_, Float_t* angle2_,
             bool printOut = false) {
    Int_t nParticles = branchParticle->GetEntries();
    Int_t levels[100];
    Int_t iLevel = 0;

    GenParticle* particle;
    GenParticle* particleM;
    // find B mesons
    for (Int_t ip = 0; ip < nParticles; ip++) {
        particle = (GenParticle*)branchParticle->At(ip);
        Int_t iMother = ip;
        Int_t haveB = 0;

        if (abs(int((particle->PID % 10000) / 100)) == 5) haveB += 1;
        if (haveB == 0) continue;

        Int_t level = 0;
        while (true) {
            particleM = (GenParticle*)branchParticle->At(iMother);
            iMother = particleM->M1;
            if (abs(particleM->PID) == 23) break;
            level += 1;
        }
        if (printOut) cout << "Found B meson: " << particle->PID << "\t(product after " << level << " decays from Z)" << endl;
        levels[iLevel] = level;
        iLevel += 1;
    }

    std::list<double> levelsList(levels, levels + iLevel);
    levelsList.sort();
    levelsList.unique();
    int levelArr[levelsList.size()];
    std::copy(levelsList.begin(), levelsList.end(), levelArr);

    Int_t idxList[999] = {0};
    Int_t nB = 0;
    Int_t nC = 0;
    if (printOut) cout << "--------------------------------------------------------------" << endl;
    for (Int_t iSize = 0; iSize < levelsList.size(); iSize++) {
        for (Int_t ip = 0; ip < nParticles; ip++) {
            string chain;
            Int_t iMother = ip;
            Int_t fromQCD = 1;
            particleM = (GenParticle*)branchParticle->At(iMother);
            for (Int_t _ = 0; _ < levelArr[iSize]; _++) {
                string chain_i = std::to_string(particleM->PID);
                iMother = particleM->M1;
                if (iMother == -1) break;
                particleM = (GenParticle*)branchParticle->At(iMother);
                chain = chain + chain_i + "\t< ";
                if (abs(int(abs(particleM->PID) / 100)) != 0) {
                    fromQCD = 0;
                    break;
                }
            }
            if (fromQCD == 0) continue;
            chain = chain + "23";

            if (iMother == -1) continue;
            particleM = (GenParticle*)branchParticle->At(iMother);
            if (abs(particleM->PID) != 23) continue;
            particle = (GenParticle*)branchParticle->At(ip);

            // c pid: 4xx, 104xx, 204xx
            if (abs(int((particle->PID % 10000) / 100)) != 5 && abs(int((particle->PID % 10000) / 100)) != 4) continue;

            if (abs(int((particle->PID % 10000) / 100)) == 5) {
                idxList[nB + nC] = ip;
                nB += 1;
            }
            if (abs(int((particle->PID % 10000) / 100)) == 4) {
                idxList[nB + nC] = ip;
                nC += 1;
            }

            TLorentzVector meson;
            meson.SetPtEtaPhiE(particle->PT, particle->Eta, particle->Phi, particle->E);
            // chain = std::to_string(level_arr[i_size]) + "-th product: (P=" + std::to_string(meson.P()) + ")\t" + chain;
            chain = std::to_string(levelArr[iSize]) + "-th product: (E=" + std::to_string(meson.E()) + ")\t" + chain;

            if (printOut) cout << chain << endl;
        }
    }
    if (printOut) cout << "--------------------------------------------------------------" << endl;
    if (printOut) cout << "nB: " << nB << "; nC: " << nC << endl;

    Float_t angles[999] = {999};
    Float_t energies[999] = {999};
    Int_t iAngles = 0;
    Int_t iEnergies = 0;
    GenParticle* particleQ1;
    GenParticle* particleQ2;
    for (int iQ1 = 0; iQ1 < nB + nC; iQ1++) {
        TLorentzVector meson1;
        particleQ1 = (GenParticle*)branchParticle->At(idxList[iQ1]);
        meson1.SetPtEtaPhiE(particleQ1->PT, particleQ1->Eta, particleQ1->Phi, particleQ1->E);
        for (int iQ2 = 0; iQ2 < iQ1; iQ2++) {
            TLorentzVector meson2;
            particleQ2 = (GenParticle*)branchParticle->At(idxList[iQ2]);
            meson2.SetPtEtaPhiE(particleQ2->PT, particleQ2->Eta, particleQ2->Phi, particleQ2->E);

            Float_t theta = acos((meson1.Px() * meson2.Px() + meson1.Py() * meson2.Py() + meson1.Pz() * meson2.Pz()) / (meson1.P() * meson2.P()));
            if (printOut) cout << "Angle Bewteen " << particleQ1->PID << " & " << particleQ2->PID << ": " << theta * 180.0 / 3.141592654 << " degrees." << endl;
            angles[iAngles] = theta * 180.0 / 3.141592654;
            iAngles += 1;
        }
        energies[iEnergies] = meson1.E();
        iEnergies += 1;
    }
    sort(angles, angles + iAngles);
    sort(energies, energies + iEnergies);

    *nB_ = nB;
    *nC_ = nC;
    *E1_ = energies[0];
    *E2_ = energies[1];
    *angle1_ = angles[0];
    *angle2_ = angles[1];
}
//}}}

//{{{
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

    Int_t nB;
    Int_t nC;
    Float_t angle1;
    Float_t angle2;
    Float_t E1;
    Float_t E2;
};
//}}}

//{{{
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
//}}}

Int_t nEvt = 0;
void eventshape(const string type, const Bool_t save = false, Int_t num_test = 0, bool printOut = false) {
    cout << "\n\n\n\n\n\n\n\n\n\n";

    string typeName;
    const char* inputFile;
    const char* outputFile;
    if (type == "s1") {
        typeName = "Bc->Jpsi tau nu. ";
        cout << "Bc->Jpsi tau nu. " << endl;
        inputFile = "./Bcjpsitaunu_50m.root";
        outputFile = "./features/JpsiTauNu_FWM.root";
    } else if (type == "s2") {
        cout << "Bc->Jpsi mu nu. " << endl;
        inputFile = "./BcJpsimunu0-2.root";
        outputFile = "./features/JpsiMuNu_FWM.root";
    } else if (type == "b1") {
        cout << "Comb+CascadeBkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        outputFile = "./features/RJpsiCombCascade_FWM.root";
    } else if (type == "b4") {
        cout << "Inclusive Bkg. " << endl;
        // inputFile = "./RJpsi_comb_200m_seed1.root";
        // outputFile = "./features/RJpsiInclusive_FWM.root";
        inputFile = "./RJpsi_comb_200m_seed2.root";
        outputFile = "./features/RJpsiInclusive_FWM_seed2.root";
    } else if (type == "b5") {
        cout << "MisID Bkg. " << endl;
        inputFile = "./RJpsi_comb_200m_seed1.root";
        outputFile = "./features/RJpsiMisID_FWM.root";
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

        if (type == "b1" || type == "b4") {
            if (printOut) cout << "\n============ Event: " << i_en << " ===================" << endl;
            Int_t nB, nC;
            Float_t E1, E2, angle1, angle2;
            calBsCs(branchParticle, &nB, &nC, &E1, &E2, &angle1, &angle2, printOut);
            FWMs->nB = nB;
            FWMs->nC = nC;
            FWMs->E1 = E1;
            FWMs->E2 = E2;
            FWMs->angle1 = angle1;
            FWMs->angle2 = angle2;
        }

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
