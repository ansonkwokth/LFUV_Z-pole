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

void mHcSpectrum(
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
        inputFile = "./Lambdab_lambdactaunu_pkpi_allchannel_1m_seed0_new.root";
        if (noise_ == 10) outputFile = "./features/mHC_Spectrum_LambdacTauNu_10Noise.root";

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
    Int_t nFoundFinalStates = 0;
    Int_t nRecoHc = 0;
    Int_t nRecoHb = 0;
    Int_t nPi_MisID = 0;  // total number of possible mis-id event (before preselection)

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

    cout << ".........." << num_test << endl;

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
        Int_t nPi_MisID_i = 0;
        // check if the event pass the classification of the correspdoning event type
        // and store the truth level info (for signal types)
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
        if (iFS.foundFromC == 0) continue;  // flag to check all of them are found

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
            nFoundFinalStates += 1;

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
            if (Length(v3Mu.X(), v3Mu.Y(), v3Mu.Z()) <= 0) continue;  // avoid the PV muon, if adding noise to make it accidently pass the cut
            v3Mu += v3MuNoise;

            // ===================================
            // ========== Lambdac cut   ==========
            // ===================================
            // c-hadron vertex distance cut
            if (Length(v3C.X(), v3C.Y(), v3C.Z()) < 0.5) continue;
            // c-hadron mass cut

            // event index
            features->iEvt = i_en;

            // pHc, EHc
            features->pHc = Hc.P();
            features->EHc = Hc.E();
            // fill the output file
            tr.Fill();

            // count the final reconstructed events
        }
    }

    if (save) {
        tr.Write();
        fea.Close();
    }
    cout << "Writing to:\t" << outputFile << endl;
}

