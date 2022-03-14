
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
//Load Delphes for reading the data
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

//---------- return if the input is signal {{{
Int_t ClassifySignal(TClonesArray* branchParticle, TLorentzVector* BTrue, TLorentzVector* HcTrue, TLorentzVector* muTrue, TVector3* v3HcTrue, TVector3* v3MuTrue, iFinalStates* iFSTrue, Int_t mode) {
    GenParticle* particle;
    GenParticle* particleB;
    GenParticle* particleC;
    GenParticle* particle0;
    GenParticle* particle0M;
    GenParticle* particle0MM;
    Int_t nParticles = branchParticle->GetEntries();

    Int_t iBTrue, foundBTrue = 0;
    Int_t iCTrue, foundCTrue = 0;
    Int_t iC1True, foundC1True = 0;
    Int_t iC2True, foundC2True = 0;
    Int_t iC3True, foundC3True = 0;
    Int_t iMuTrue, foundMuTrue = 0;
    for (int ipB = 0; ipB < nParticles; ipB++) {
        foundBTrue = 0;
        particleB = (GenParticle*)branchParticle->At(ipB);
        if (abs(particleB->PID) == 531) {
            iBTrue = ipB;
            foundBTrue = 1;

            Int_t iCTrue_i, foundCTrue_i = 0;    // Ds
            Int_t iC0True_i, foundC0True_i = 0;  // phi
            Int_t iC1True_i, foundC1True_i = 0;  // K-
            Int_t iC2True_i, foundC2True_i = 0;  // K+
            Int_t iC3True_i, foundC3True_i = 0;  // pi
            Int_t iMuTrue_i, foundMuTrue_i = 0;  // mu
            for (int ip2 = 0; ip2 < nParticles; ip2++) {
                particle0 = (GenParticle*)branchParticle->At(ip2);
                if (particle0->M1 == -1) {
                    continue;
                }
                particle0M = (GenParticle*)branchParticle->At(particle0->M1);
                if (particle0M->M1 == -1) {
                    continue;
                }
                particle0MM = (GenParticle*)branchParticle->At(particle0M->M1);

                if (abs(particle0->PID) == 431 && particle0->M1 == iBTrue) {
                    iCTrue_i = ip2;
                    foundCTrue_i = 1;
                }
                if (abs(particle0->PID) == 333 && abs(particle0M->PID) == 431 &&
                    particle0M->M1 == iBTrue) {
                    iC0True_i = ip2;
                    foundC0True_i = 1;
                }
                if (particle0->PID == -321 && abs(particle0M->PID) == 333 &&
                    abs(particle0MM->PID) == 431 && particle0MM->M1 == iBTrue) {
                    iC1True_i = ip2;
                    foundC1True_i = 1;
                }
                if (particle0->PID == 321 && abs(particle0M->PID) == 333 &&
                    abs(particle0MM->PID) == 431 && particle0MM->M1 == iBTrue) {
                    iC2True_i = ip2;
                    foundC2True_i = 1;
                }
                if (abs(particle0->PID) == 211 && abs(particle0M->PID) == 431 &&
                    particle0M->M1 == iBTrue) {
                    iC3True_i = ip2;
                    foundC3True_i = 1;
                }

                if (mode == 1 &&
                    (abs(particle0->PID) == 13 && abs(particle0M->PID) == 15 &&
                     particle0M->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
                if (mode == 2 &&
                    (abs(particle0->PID) == 13 && particle0->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
            }
            // cout <<  foundCTrue_i <<"; "<< foundC1True_i <<"; "<<
            // foundC2True_i <<"; "<< foundC3True_i<< "; " <<
            // foundMuTrue_i<<endl;
            if (foundCTrue_i == 1 && foundC1True_i == 1 && foundC2True_i == 1 &&
                foundC3True_i == 1 && foundMuTrue_i == 1) {
                iCTrue = iCTrue_i;
                foundCTrue = 1;
                iC1True = iC1True_i;
                foundC1True = 1;
                iC2True = iC2True_i;
                foundC2True = 1;
                iC3True = iC3True_i;
                foundC3True = 1;
                iMuTrue = iMuTrue_i;
                foundMuTrue = 1;
                break;
            }
        }
    }
    // cout <<  foundBTrue <<"; "<< foundCTrue <<"; "<< foundCTrue <<"; "<<
    // foundC2True <<"; "<< foundC3True << "; " << foundMuTrue <<endl;
    if (foundBTrue == 1 && foundCTrue == 1 && foundC1True == 1 &&
        foundC2True == 1 && foundC3True == 1 && foundMuTrue == 1) {
        iFinalStates iFSTrue_;
        iFSTrue_.iKNeg = iC1True;
        iFSTrue_.iKPos = iC2True;
        iFSTrue_.iPi = iC3True;
        iFSTrue_.iMu = iMuTrue;

        TLorentzVector BTrue_, HcTrue_, muTrue_;
        GenParticle* particleHc;
        GenParticle* particleMu;
        GenParticle* particleC1;
        particleB = (GenParticle*)branchParticle->At(iBTrue);
        BTrue_.SetPtEtaPhiE(particleB->PT, particleB->Eta, particleB->Phi, particleB->E);
        particleHc = (GenParticle*)branchParticle->At(iCTrue);
        HcTrue_.SetPtEtaPhiE(particleHc->PT, particleHc->Eta, particleHc->Phi, particleHc->E);
        particleMu = (GenParticle*)branchParticle->At(iMuTrue);
        muTrue_.SetPtEtaPhiE(particleMu->PT, particleMu->Eta, particleMu->Phi, particleMu->E);

        particleC1 = (GenParticle*)branchParticle->At(iC1True);
        TVector3 v3HcTrue_, v3MuTrue_;
        v3HcTrue_.SetXYZ(particleC1->X, particleC1->Y, particleC1->Z);
        v3MuTrue_.SetXYZ(particleMu->X, particleMu->Y, particleMu->Z);

        *iFSTrue = iFSTrue_;
        *BTrue = BTrue_;
        *HcTrue = HcTrue_;
        *muTrue = muTrue_;
        *v3HcTrue = v3HcTrue_;
        *v3MuTrue = v3MuTrue_;
        return 1;
    } else {
        return 0;
    }
}
//}}}

//---------- return if the input is signal {{{
Int_t ClassifyExcitedSignal(TClonesArray* branchParticle, TLorentzVector* BTrue, TLorentzVector* HcTrue, TLorentzVector* muTrue, TLorentzVector* phoTrue, TVector3* v3HcTrue, TVector3* v3MuTrue, iFinalStates* iFSTrue, Int_t mode) {
    GenParticle* particle;
    GenParticle* particleB;
    GenParticle* particleC;
    GenParticle* particlePho;
    GenParticle* particle0;
    GenParticle* particle0M;
    GenParticle* particle0MM;
    GenParticle* particle0MMM;
    Int_t nParticles = branchParticle->GetEntries();

    Int_t iBTrue, foundBTrue = 0;
    Int_t iCTrue, foundCTrue = 0;
    Int_t iPhoTrue, foundPhoTrue = 0;
    Int_t iC1True, foundC1True = 0;
    Int_t iC2True, foundC2True = 0;
    Int_t iC3True, foundC3True = 0;
    Int_t iMuTrue, foundMuTrue = 0;
    for (int ipB = 0; ipB < nParticles; ipB++) {
        foundBTrue = 0;
        particleB = (GenParticle*)branchParticle->At(ipB);
        if (abs(particleB->PID) == 531) {
            iBTrue = ipB;
            foundBTrue = 1;

            Int_t iCTrue_i, foundCTrue_i = 0;    // Ds
            Int_t iC0True_i, foundC0True_i = 0;  // phi
            Int_t iC1True_i, foundC1True_i = 0;  // K-
            Int_t iC2True_i, foundC2True_i = 0;  // K+
            Int_t iC3True_i, foundC3True_i = 0;  // pi
            Int_t iMuTrue_i, foundMuTrue_i = 0;  // mu
            for (int ip2 = 0; ip2 < nParticles; ip2++) {
                particle0 = (GenParticle*)branchParticle->At(ip2);
                if (particle0->M1 == -1) {
                    continue;
                }
                particle0M = (GenParticle*)branchParticle->At(particle0->M1);
                if (particle0M->M1 == -1) {
                    continue;
                }
                particle0MM = (GenParticle*)branchParticle->At(particle0M->M1);
                if (particle0MM->M1 == -1) {
                    continue;
                }
                particle0MMM =
                    (GenParticle*)branchParticle->At(particle0MM->M1);

                if (abs(particle0->PID) == 431 && abs(particle0M->PID) == 433 &&
                    particle0M->M1 == iBTrue) {
                    iCTrue_i = ip2;
                    foundCTrue_i = 1;
                }
                if (abs(particle0->PID) == 333 && abs(particle0->PID) == 431 &&
                    abs(particle0M->PID) == 433 && particle0MM->M1 == iBTrue) {
                    iC0True_i = ip2;
                    foundC0True_i = 1;
                }
                if (particle0->PID == -321 && abs(particle0M->PID) == 333 &&
                    abs(particle0MM->PID) == 431 &&
                    abs(particle0MMM->PID) == 433 &&
                    particle0MMM->M1 == iBTrue) {
                    iC1True_i = ip2;
                    foundC1True_i = 1;
                }
                if (particle0->PID == 321 && abs(particle0M->PID) == 333 &&
                    abs(particle0MM->PID) == 431 &&
                    abs(particle0MMM->PID) == 433 &&
                    particle0MMM->M1 == iBTrue) {
                    iC2True_i = ip2;
                    foundC2True_i = 1;
                }
                if (abs(particle0->PID) == 211 && abs(particle0M->PID) == 431 &&
                    abs(particle0MM->PID) == 433 && particle0MM->M1 == iBTrue) {
                    iC3True_i = ip2;
                    foundC3True_i = 1;
                }

                if (mode == 1 &&
                    (abs(particle0->PID) == 13 && abs(particle0M->PID) == 15 &&
                     particle0M->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
                if (mode == 2 &&
                    (abs(particle0->PID) == 13 && particle0->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
            }
            // cout <<  foundCTrue_i <<"; "<< foundC1True_i <<"; "<<
            // foundC2True_i <<"; "<< foundC3True_i<< "; " <<
            // foundMuTrue_i<<endl;
            if (foundCTrue_i == 1 && foundC1True_i == 1 && foundC2True_i == 1 &&
                foundC3True_i == 1 && foundMuTrue_i == 1) {
                iCTrue = iCTrue_i;
                foundCTrue = 1;
                iC1True = iC1True_i;
                foundC1True = 1;
                iC2True = iC2True_i;
                foundC2True = 1;
                iC3True = iC3True_i;
                foundC3True = 1;
                iMuTrue = iMuTrue_i;
                foundMuTrue = 1;
                break;
            }
        }
    }

    if (not(foundBTrue == 1 && foundCTrue == 1 && foundCTrue == 1 &&
            foundC2True == 1 && foundC3True == 1 && foundMuTrue == 1)) {
        return 0;
    }
    for (int ipho = 0; ipho < nParticles; ipho++) {
        particlePho = (GenParticle*)branchParticle->At(ipho);
        particleC = (GenParticle*)branchParticle->At(iCTrue);
        if (abs(particlePho->PID) != 22) {
            continue;
        }
        if (particlePho->M1 == particleC->M1) {
            foundPhoTrue = 1;
            iPhoTrue = ipho;
        }
    }

    // cout<< foundPhoTrue << "; "<<  foundBTrue <<"; "<< foundCTrue <<"; "<<
    // foundCTrue <<"; "<< foundC2True <<"; "<< foundC3True << "; " <<
    // foundMuTrue <<endl;
    if (foundPhoTrue == 1 && foundBTrue == 1 && foundCTrue == 1 &&
        foundCTrue == 1 && foundC2True == 1 && foundC3True == 1 &&
        foundMuTrue == 1) {
        iFinalStates iFSTrue_;
        iFSTrue_.iKNeg = iC1True;
        iFSTrue_.iKPos = iC2True;
        iFSTrue_.iPi = iC3True;
        iFSTrue_.iMu = iMuTrue;

        TLorentzVector BTrue_, HcTrue_, muTrue_, phoTrue_;
        GenParticle* particleHc;
        GenParticle* particleMu;
        GenParticle* particleC1;
        particleB = (GenParticle*)branchParticle->At(iBTrue);
        BTrue_.SetPtEtaPhiE(particleB->PT, particleB->Eta, particleB->Phi, particleB->E);
        particleHc = (GenParticle*)branchParticle->At(iCTrue);
        HcTrue_.SetPtEtaPhiE(particleHc->PT, particleHc->Eta, particleHc->Phi, particleHc->E);
        particleMu = (GenParticle*)branchParticle->At(iMuTrue);
        muTrue_.SetPtEtaPhiE(particleMu->PT, particleMu->Eta, particleMu->Phi, particleMu->E);
        particlePho = (GenParticle*)branchParticle->At(iPhoTrue);
        phoTrue_.SetPtEtaPhiE(particlePho->PT, particlePho->Eta, particlePho->Phi, particlePho->E);

        particleC1 = (GenParticle*)branchParticle->At(iC1True);
        TVector3 v3HcTrue_, v3MuTrue_;
        v3HcTrue_.SetXYZ(particleC1->X, particleC1->Y, particleC1->Z);
        v3MuTrue_.SetXYZ(particleMu->X, particleMu->Y, particleMu->Z);

        *iFSTrue = iFSTrue_;
        *BTrue = BTrue_;
        *HcTrue = HcTrue_;
        *muTrue = muTrue_;
        *phoTrue = phoTrue_;
        *v3HcTrue = v3HcTrue_;
        *v3MuTrue = v3MuTrue_;
        return 1;
    } else {
        return 0;
    }
}
//}}}

//---------- return if the event is the targeted bkg (from inclusive samples)
//{{{
Int_t ClassifyBkg(TClonesArray* branchParticle, const string type) {
    // Inclusive samepls (for Comb+Cascade, Comb, Casde, Inclusive)
    GenParticle* particle;
    Int_t nParticles = branchParticle->GetEntries();

    //=== check c-hadron ===
    Int_t nPosHc = 0;  // number of positive targeted c-hadron in the truth
                       // level (Jpsi, Ds, Lambdac)
    Int_t nNegHc = 0;  // number of negative targeted c-hadron in the truth
                       // level (Jpsi, Ds, Lambdac)
    for (int ip = 0; ip < nParticles; ip++) {
        particle = (GenParticle*)branchParticle->At(ip);
        if (particle->PID == 431) {
            nPosHc += 1;
        }  // check number of pos lambdac
        if (particle->PID == -431) {
            nNegHc += 1;
        }  // check number of neg lambdac
    }
    if (not(nPosHc + nNegHc >= 1)) {
        return 0;
    }  // must have at least one lambdac

    //=== check muons ===
    Int_t nPosMu = 0;  // number of positive muons
    Int_t nNegMu = 0;  // number of negative muons
    for (int ip = 0; ip < nParticles; ip++) {
        particle = (GenParticle*)branchParticle->At(ip);
        if (particle->PID == -13) {
            nPosMu += 1;
        }
        if (particle->PID == 13) {
            nNegMu += 1;
        }
    }
    // total must be at least 3, (2 from Jpsi, 1 unapired), and at least 1 pos,
    // 1 neg from Jpsi
    if (not((nPosHc > 0 && nNegMu > 0) || (nNegHc > 0 && nPosMu > 0))) {
        return 0;
    }

    //=== check the decay chain ===
    Int_t isLeptonFromB = 0;  // check if the lepton is from b-hadron
    Int_t isLeptonSemiFromB =
        0;  // check if the lepton is semileptonically decay from b-hadron
    Int_t nLeptonStepFromB = 0;
    Int_t isCFromB = 0;  // check if c-hadron is from B
    Int_t isSignal = 0;  // chekc if event has signal in it

    GenParticle* particle1;
    GenParticle* particle1M;  // particle1's mother
    GenParticle* particle2;
    GenParticle* particle2M;   // particle2's mother
    GenParticle* particle2MM;  // particle2's mother
    GenParticle* particleL;    // lepton
    GenParticle* particleC;    // c-hadron
    for (int ip1 = 0; ip1 < nParticles; ip1++) {
        particle1 = (GenParticle*)branchParticle->At(ip1);
        // identify the muon
        if (abs(particle1->PID) == 13) {
            isLeptonFromB = 0;
            isLeptonSemiFromB = 0;
            isCFromB = 0;

            // mother of particle0
            particle1M = (GenParticle*)branchParticle->At(particle1->M1);

            // identify the lepton (mu/tau)
            Int_t iLepton = ip1;
            if (abs(particle1M->PID) == 15) {
                iLepton = particle1->M1;
            }
            particleL = (GenParticle*)branchParticle->At(iLepton);
            Int_t iLeptonMother = particleL->M1;
            particle1M = (GenParticle*)branchParticle->At(iLeptonMother);
            // cout << "Lepton: " << particleL->PID << " (" << iLepton << ")
            // from: ";

            // finding the b-hadron that lepton from
            Int_t BHadron_idx = 99999;
            nLeptonStepFromB = 0;
            while (true) {
                // cout << particle1M->PID << " < ";
                if (int(abs(particle1M->PID) / 100) == 5 ||
                    int(abs(particle1M->PID) / 1000) == 5 ||
                    int(abs(particle1M->PID) / 10000) == 5 ||
                    int((abs(particle1M->PID) % 1000) / 100) == 5) {
                    BHadron_idx = iLeptonMother;
                    // cout << "\n From b-hadron (IDX): " << particle1M->PID <<
                    // " (" << BHadron_idx << ")" << endl;
                    break;
                } else if (int(particle1M->PID / 10) == 0) {
                    // cout << "\nNo b-hadron, terminated at: " <<
                    // particle1M->PID << endl;
                    break;
                }
                nLeptonStepFromB += 1;
                iLeptonMother = particle1M->M1;
                particle1M = (GenParticle*)branchParticle->At(iLeptonMother);
            }
            if (BHadron_idx != 99999) {
                isLeptonFromB = 1;
            }
            if (BHadron_idx != 99999 && nLeptonStepFromB == 0) {
                isLeptonSemiFromB = 1;
            }

            // if lepton is from b-hadron, then find if c-hadron is from the
            // same mother
            if (isLeptonFromB == 1) {
                for (int ip2 = 0; ip2 < nParticles; ip2++) {
                    particle2 = (GenParticle*)branchParticle->At(ip2);
                    if (abs(particle2->PID) == 431) {
                        Int_t iC = ip2;
                        particleC = (GenParticle*)branchParticle->At(iC);
                        Int_t iCMother = particleC->M1;
                        particle2M =
                            (GenParticle*)branchParticle->At(iCMother);
                        Int_t iCMotherMother = particle2M->M1;
                        if (particle2M->M1 == -1) {
                            continue;
                        }
                        particle2MM =
                            (GenParticle*)branchParticle->At(particle2M->M1);

                        if (iCMother == BHadron_idx && isLeptonSemiFromB == 1) {
                            isSignal = 1;
                            continue;
                        }
                        if (abs(particle2M->PID) == 433 &&
                            iCMotherMother == BHadron_idx &&
                            isLeptonSemiFromB == 1) {
                            isSignal = 1;
                            continue;
                        }
                        // cout << "c-hadron " << particleC->PID << " (" <<
                        // iCMother << ") from: ";

                        while (true) {
                            // cout << particle2M->PID << " (" << iCMother <<
                            // ")" << " < ";
                            if (iCMother == BHadron_idx) {
                                // cout << "/nc-hadron from the smae b-hadron: "
                                // << particle2M->PID << endl;
                                isCFromB = 1;
                                break;
                            }
                            if (int(particle2M->PID / 10) == 0) {
                                // cout << "\nc-hadron not from same b-hadron
                                // terminate at: " << particle2M->PID << endl;
                                break;
                            }
                            iCMother = particle2M->M1;
                            particle2M =
                                (GenParticle*)branchParticle->At(iCMother);
                        }
                        if (isCFromB == 1 || isSignal == 1) {
                            break;
                        }
                    }
                }
                if (isCFromB == 1 || isSignal == 1) {
                    break;
                }
            }
        }
    }
    if (isSignal == 1) {
        return 0;
    }

    Int_t isComb = 0;
    Int_t isCascade = 0;
    Int_t isInclusive = 0;

    if (isLeptonSemiFromB == 1 && isCFromB == 1) {
        isInclusive = 1;
    }
    if (isLeptonFromB == 1 && isLeptonSemiFromB == 0 && isCFromB == 1) {
        isCascade = 1;
    }
    if (isCFromB == 0) {
        isComb = 1;
    }

    if (isInclusive + isCascade + isComb != 1) {
        cout << "HAVE BUG IN CLASSIFYING BKG!";
    }

    // cout << "Signal: " << isSignal << "; Inclsuive: " << isInclusive << ";
    // Cascade: " << isCascade << "; Comb:" << isComb << endl; cout << endl;
    if (type == "b1") {
        if (isCascade == 0 && isComb == 0) {
            return 0;
        } else {
            return 1;
        }
    }
    if (type == "b2") {
        if (isComb == 0) {
            return 0;
        } else {
            return 1;
        }
    }
    if (type == "b3") {
        if (isCascade == 0) {
            return 0;
        } else {
            return 1;
        }
    }
    if (type == "b4") {
        if (isInclusive == 0) {
            return 0;
        } else {
            return 1;
        }
    } else {
        return 0;
    }
}
//}}}

//---------- return if the event is misID bkg{{{
Int_t ClassifyMisID(TClonesArray* branchParticle) {
    // skip signal events
    TLorentzVector dummy;
    TVector3 dummy2;
    iFinalStates dummy3;
    if (ClassifySignal(branchParticle, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy3, 1) == 1) {
        return 0;
    } else if (ClassifySignal(branchParticle, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy3, 2) == 1) {
        return 0;
    } else if (ClassifyExcitedSignal(branchParticle, &dummy, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy3, 1) == 1) {
        return 0;
    } else if (ClassifyExcitedSignal(branchParticle, &dummy, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy3, 2) == 1) {
        return 0;
    } else {
        // cout<< 1 <<endl;
        GenParticle* particle;
        Int_t nParticles = branchParticle->GetEntries();

        //=== check c-hadron ===
        Int_t nHc = 0;  // number of targeted c-hadron in the truth level (Jpsi,
                        // Ds, Lambdac)
        for (int ip = 0; ip < nParticles; ip++) {
            particle = (GenParticle*)branchParticle->At(ip);
            if (abs(particle->PID) == 431) {
                nHc += 1;
            }  // check number of Jpsi
        }
        if (not(nHc >= 1)) {
            return 0;
        }  // must have at least one Jpsi

        //=== check pions ===
        Int_t nPi = 0;  // number of pions (displaced)
        for (int ip = 0; ip < nParticles; ip++) {
            particle = (GenParticle*)branchParticle->At(ip);
            if (abs(particle->PID) == 211 &&
                (not(particle->X == 0 && particle->Y == 0 &&
                     particle->Z == 0))) {
                nPi += 1;
            }
        }
        if (not(nPi >= 2)) {
            return 0;
        }
        return 1;
    }
}
//}}}

void allinone(const string type, const Float_t noise_ = 10, const Bool_t save = false, Int_t num_test = 0, const Float_t alpha = 1) {
    cout << "\n\n\n\n\n\n\n\n"
         << endl;
    string typeName;
    const char* inputFile;
    const char* outputFile;
    const Float_t noise = noise_ * 0.001;
    if (type == "s1") {
        typeName = "Bs->Ds tau nu. ";
        cout << "Bs->Ds tau nu. " << endl;
        inputFile = "./Bs0Dstaunu-Dsphipi-phiKK_100k_RandomSeed0.root";
        if (noise_ == 10) {
            outputFile = "./features/DsTauNu_10Noise.root";
        } else if (noise_ == 20) {
            outputFile = "./features/DsTauNu_20Noise.root";
        }

    } else if (type == "s2") {
        cout << "Bs->Ds mu nu. " << endl;
        inputFile = "./Bs0Dsmunu-Dsphipi-phiKK_100k_RandomSeed0.root";
        if (noise_ == 10) {
            outputFile = "./features/DsMuNu_10Noise.root";
        } else if (noise_ == 20) {
            outputFile = "./features/DsMuNu_20Noise.root";
        }
    } else if (type == "s3") {
        cout << "Bs->Ds* mu nu. " << endl;
        inputFile = "./Bs0Dsstartaunu-Dsphipi-phiKK_30k_RandomSeed0_30k_RandomSeed1.root";
        if (noise_ == 10) {
            outputFile = "./features/DsstarTauNu_10Noise.root";
            if (abs(alpha - 0) < 1e-6) {
                outputFile = "./features/DsstarTauNu_10Noise_0Alpha.root";
            } else if (abs(alpha - 0.5) < 1e-6) {
                outputFile = "./features/DsstarTauNu_10Noise_05Alpha.root";
            } else if (abs(alpha - 2) < 1e-6) {
                outputFile = "./features/DsstarTauNu_10Noise_2Alpha.root";
            } else if (abs(alpha - 0.1) < 1e-6) {
                outputFile = "./features/DsstarTauNu_10Noise_01Alpha.root";
            }
        } else if (noise_ == 20) {
            outputFile = "./features/DsstarTauNu_20Noise.root";
        }
    } else if (type == "s4") {
        cout << "Bs->Ds* mu nu. " << endl;
        inputFile = "./Bs0Dsstarmunu-Dsphipi-phiKK_30k_RandomSeed0.root";
        if (noise_ == 10) {
            outputFile = "./features/DsstarMuNu_10Noise.root";
            if (abs(alpha - 0) < 1e-6) {
                outputFile = "./features/DsstarMuNu_10Noise_0Alpha.root";
            } else if (abs(alpha - 0.1) < 1e-6) {
                outputFile = "./features/DsstarMuNu_10Noise_01Alpha.root";
            } else if (abs(alpha - 0.5) < 1e-6) {
                outputFile = "./features/DsstarMuNu_10Noise_05Alpha.root";
            } else if (abs(alpha - 2) < 1e-6) {
                outputFile = "./features/DsstarMuNu_10Noise_2Alpha.root";
            }
        } else if (noise_ == 20) {
            outputFile = "./features/DsstarMuNu_20Noise.root";
        }
    } else if (type == "b1") {
        cout << "Comb+Cascade Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        outputFile = "./testoutput";
        if (noise_ == 10) {
            outputFile = "./features/RDsCombCascade_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/RDsCombCascade_20Noise.root";
        }
    } else if (type == "b2") {
        cout << "Comb Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        outputFile = "./testoutput";
        if (noise_ == 10) {
            outputFile = "./features/RDsComb_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/RDsComb_20Noise.root";
        }
    } else if (type == "b3") {
        cout << "Cascade Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        outputFile = "./testoutput";
        if (noise_ == 10) {
            outputFile = "./features/RDsCascade_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/RDsCascade_20Noise.root";
        }
    } else if (type == "b4") {
        cout << "Inclusive Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        outputFile = "./testoutput.root";
        if (noise_ == 10) {
            outputFile = "./features/RDsInclusive_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/RDsInclusive_20Noise.root";
        }
    } else if (type == "b5") {
        cout << "MisID Bkg. " << endl;
        inputFile = "./RDs_comb_1mseed1_1mseed2_1mseed3.root";
        outputFile = "./testoutput";
        if (noise_ == 10) {
            outputFile = "./features/RDsMisID_10Noise.root";
        }
        if (noise_ == 20) {
            outputFile = "./features/RDsMisID_20Noise.root";
        }
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
    TClonesArray* branchEFlowNeutralHadron =
        treeReader->UseBranch("EFlowNeutralHadron");

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

    GenParticle* particle;
    Track* track;
    Tower* eflowphoton;
    Tower* eflowneutralhadron;

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
    Features* features = new Features;
    // Features features;
    tr.Branch("features", &features);

    cout << "\nReconstructing: " << typeName << endl;
    cout << "Reading from:\t" << inputFile << endl;
    cout << "Injected noise:\t" << noise_ << " microns (" << noise << "mm)"
         << endl;
    cout << "Save: " << save << endl;
    cout << "Writing to:\t" << outputFile << endl;
    cout << "# of events:\t" << numberOfEntries << "\n"
         << endl;
    if (num_test == 0) {
        num_test = numberOfEntries;
    }

    // for (Int_t i_en=0; i_en< numberOfEntries; i_en++) {
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // if (i_en % 1000 == 0) {cout << "\tReconstruction Progress: " << i_en
        // << "/" << numberOfEntries; } if ((i_en % 1000) == 0) {cout <<
        // "\rReconstruction Progress: " << i_en << "/" << numberOfEntries; }
        if ((i_en % 100000) == 0) {
            cout << "Reconstruction Progress: " << i_en << "/"
                 << numberOfEntries << endl;
        }
        // if (i_en >= num_test) { break; }
        treeReader->ReadEntry(i_en);  // reading the entry

        //==========================================================================
        //===============   Classifying in event type in truth level
        //===============
        //==========================================================================
        iFinalStates iFSTrue;
        TLorentzVector BTrue, HcTrue, muTrue, phoTrue;
        TVector3 v3HcTrue, v3MuTrue;
        Int_t passing = 0;
        if (type == "s1") {
            passing = ClassifySignal(branchParticle, &BTrue, &HcTrue, &muTrue, &v3HcTrue, &v3MuTrue, &iFSTrue, 1);
        } else if (type == "s2") {
            passing = ClassifySignal(branchParticle, &BTrue, &HcTrue, &muTrue, &v3HcTrue, &v3MuTrue, &iFSTrue, 2);
        } else if (type == "s3") {
            passing = ClassifyExcitedSignal(branchParticle, &BTrue, &HcTrue, &muTrue, &phoTrue, &v3HcTrue, &v3MuTrue, &iFSTrue, 1);
        } else if (type == "s4") {
            passing = ClassifyExcitedSignal(branchParticle, &BTrue, &HcTrue, &muTrue, &phoTrue, &v3HcTrue, &v3MuTrue, &iFSTrue, 2);
        } else if (type == "b1" || type == "b1" || type == "b2" ||
                   type == "b3" || type == "b4") {
            passing = ClassifyBkg(branchParticle, type);
        } else if (type == "b5") {
            passing = ClassifyMisID(branchParticle);
        } else {
            cout << " No such Signal/Bkg" << endl;
        }

        // cout<<i_en <<endl;
        if (passing == 0) {
            continue;
        }
        // cout << "\n===========================================Event: " <<
        // i_en << endl;
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

        Int_t nRecoMu = 0;
        Int_t nMuCand = 0;
        for (Int_t imu : muCandidates) {
            if (imu == 99999) {
                continue;
            }
            // cout << iloop << endl;
            iFS.iMu = imu;
            // cout << " ..." << iFS.iKPos << "; " << iFS.iKNeg << "; " <<
            // iFS.iPi << "; " << iFS.iMu << endl;

            //==========================================================================
            //==================   Defining 4 Momentum & 3 Positions  ==================
            //==========================================================================

            Track* KNegTrack = (Track*)branchTrack->At(iFS.iKNeg);
            Track* KPosTrack = (Track*)branchTrack->At(iFS.iKPos);
            Track* piTrack = (Track*)branchTrack->At(iFS.iPi);
            Track* muTrack = (Track*)branchTrack->At(iFS.iMu);

            TLorentzVector KNeg;
            KNeg.SetPtEtaPhiM(KNegTrack->PT, KNegTrack->Eta, KNegTrack->Phi, mK);
            TLorentzVector KPos;
            KPos.SetPtEtaPhiM(KPosTrack->PT, KPosTrack->Eta, KPosTrack->Phi, mK);
            TLorentzVector pi;
            pi.SetPtEtaPhiM(piTrack->PT, piTrack->Eta, piTrack->Phi, mPi);
            TLorentzVector mu;
            mu.SetPtEtaPhiM(muTrack->PT, muTrack->Eta, muTrack->Phi, mMu);

            if (KNeg.Px() * KPos.Px() + KNeg.Py() * KPos.Py() + KNeg.Pz() * KPos.Pz() <= 0) {
                continue;
            }
            if (KNeg.Px() * KNeg.Px() + KNeg.Py() * KNeg.Py() + KNeg.Pz() * KNeg.Pz() <= 0) {
                continue;
            }
            if (KNeg.Px() * mu.Px() + KNeg.Py() * mu.Py() + KNeg.Pz() * mu.Pz() <= 0) {
                continue;
            }
            if (KPos.Px() * pi.Px() + KPos.Py() * pi.Py() + KPos.Pz() * pi.Pz() <= 0) {
                continue;
            }
            if (KPos.Px() * mu.Px() + KPos.Py() * mu.Py() + KPos.Pz() * mu.Pz() <= 0) {
                continue;
            }
            if (pi.Px() * mu.Px() + pi.Py() * mu.Py() + pi.Pz() * mu.Pz() <= 0) {
                continue;
            }

            nSameDir += 1;

            TLorentzVector Hc;
            Hc = KNeg + KPos + pi;

            Float_t dxC = distribution(genertator);
            Float_t dyC = distribution(genertator);
            Float_t dzC = distribution(genertator);
            Float_t dxMu = distribution(genertator);
            Float_t dyMu = distribution(genertator);
            Float_t dzMu = distribution(genertator);
            TVector3 v3CNoise(dxC, dyC, dzC);
            TVector3 v3MuNoise(dxMu, dyMu, dzMu);
            TVector3 v3C(piTrack->X, piTrack->Y, piTrack->Z);
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
            //=============================   Apply cuts
            //=============================
            //==========================================================================
            if (not(Hc.M() >= 1.945 && Hc.M() <= 1.995)) {
                continue;
            }
            nHcMass += 1;

            TLorentzVector KK;
            KK = KPos + KNeg;
            if (not(KK.M() >= 1.008 && KK.M() <= 1.032)) {
                continue;
            }

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

            distance_2lines(v3C.X(), v3C.Y(), v3C.Z(), Hc.Px(), Hc.Py(), Hc.Pz(), v3Mu.X(), v3Mu.Y(), v3Mu.Z(), mu.Px(), mu.Py(), mu.Pz(), &LHcMu, &sHc, &sMu);
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
            B = reconstructB4Momentum(v3B, iFS, branchTrack, branchEFlowPhoton, branchEFlowNeutralHadron);
            if (isnan(B.P())) {
                continue;
            }
            if (B.Px() == 99999 && B.Py() == 99999 && B.Pz() == 99999 &&
                B.E() == 99999) {
                continue;
            }

            TLorentzVector HcMu;
            HcMu = KNeg + KPos + pi + mu;
            if (HcMu.M() > mBs) {
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

            TLorentzVector qTrue;
            qTrue = BTrue - HcTrue;
            TLorentzVector missTrue;
            missTrue = BTrue - HcTrue - muTrue;
            Float_t q2True = qTrue.E() * qTrue.E() - qTrue.P() * qTrue.P();
            Float_t miss2True =
                missTrue.E() * missTrue.E() - missTrue.P() * missTrue.P();

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
            // cout << "phoE: " << features->EPhoTrue << "; " << features->correctPhoton << endl;

            tr.Fill();
            num += 1;
            nRecoMu += 1;
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

