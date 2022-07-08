//
// Testing the Tower branch and EFlowPhoton branch relation
//
//  ==================
//  === R(D_s^(*)) ===
//  ==================
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

void testTower(
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
        if (noise_ == 0) outputFile = "./features/DsTauNu_0Noise.root";
        if (noise_ == 5) outputFile = "./features/DsTauNu_5Noise.root";
        if (noise_ == 10) outputFile = "./features/DsTauNu_10Noise.root";
        if (noise_ == 20) outputFile = "./features/DsTauNu_20Noise.root";

    } else if (type == "s2") {
        cout << "Bs->Ds mu nu. " << endl;
        inputFile = "./Bs0Dsmunu-Dsphipi-phiKK_100k_RandomSeed0.root";
        if (noise_ == 0) outputFile = "./features/DsMuNu_0Noise.root";
        if (noise_ == 5) outputFile = "./features/DsMuNu_5Noise.root";
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
        } else if (noise_ == 5) {
            outputFile = "./features/DsstarTauNu_5Noise.root";
        } else if (noise_ == 0) {
            outputFile = "./features/DsstarTauNu_0Noise.root";
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
        } else if (noise_ == 5) {
            outputFile = "./features/DsstarMuNu_5Noise.root";
        } else if (noise_ == 0) {
            outputFile = "./features/DsstarMuNu_0Noise.root";
        }

    } else if (type == "b1") {
        cout << "Comb+Cascade Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        if (noise_ == 0) outputFile = "./features/RDsCombCascade_0Noise.root";
        if (noise_ == 5) outputFile = "./features/RDsCombCascade_5Noise.root";
        if (noise_ == 10) outputFile = "./features/RDsCombCascade_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RDsCombCascade_20Noise.root";

    } else if (type == "b2") {
        cout << "Comb Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        if (noise_ == 0) outputFile = "./features/RDsComb_0Noise.root";
        if (noise_ == 5) outputFile = "./features/RDsComb_5Noise.root";
        if (noise_ == 10) outputFile = "./features/RDsComb_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RDsComb_20Noise.root";
        // if (noise_ == 0) outputFile = "./features/RDsComb_0Noise.root";

    } else if (type == "b3") {
        cout << "Cascade Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        if (noise_ == 0) outputFile = "./features/RDsCascade_0Noise.root";
        if (noise_ == 5) outputFile = "./features/RDsCascade_5Noise.root";
        if (noise_ == 10) outputFile = "./features/RDsCascade_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RDsCascade_20Noise.root";
        // if (noise_ == 0) outputFile = "./features/RDsCascade_0Noise.root";

    } else if (type == "b4") {
        cout << "Inclusive Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        if (noise_ == 0) outputFile = "./features/RDsInclusive_0Noise.root";
        if (noise_ == 5) outputFile = "./features/RDsInclusive_5Noise.root";
        if (noise_ == 10) outputFile = "./features/RDsInclusive_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RDsInclusive_20Noise.root";

    } else if (type == "b5") {
        cout << "MisID Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        if (noise_ == 0) outputFile = "./features/RDsMisID_0Noise.root";
        if (noise_ == 5) outputFile = "./features/RDsMisID_5Noise.root";
        if (noise_ == 10) outputFile = "./features/RDsMisID_10Noise.root";
        if (noise_ == 20) outputFile = "./features/RDsMisID_20Noise.root";

    } else {
        cout << "Not match. ";
        return;
    }

    //
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
    TClonesArray* branchTower = treeReader->UseBranch("Tower");

    GenParticle* particle;
    Track* track;
    Tower* eflowphoton;
    Tower* tower;
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

        Int_t nEF = branchEFlowPhoton->GetEntries();
        Int_t nT = branchTower->GetEntries();
        for (Int_t ief = 0; ief < nEF; ief++) {
            eflowphoton = (Tower*)branchEFlowPhoton->At(ief);
            cout << " E: " << eflowphoton->E << "\n";
            cout << " UID: " << eflowphoton->GetUniqueID() << "\n";
            // for (Int_t it = 0; it < nT; it++) {
            // tower = (Tower*)branchTower->At(it);
            // cout << " tower UID: " << tower->GetUniqueID() << "\n";
            // if (tower->GetUniqueID() == eflowphoton->GetUniqueID()) {
            // cout << " tower E: " << tower->E << "\n";
            //}
            //}
        }
        cout << endl;
    }
}

