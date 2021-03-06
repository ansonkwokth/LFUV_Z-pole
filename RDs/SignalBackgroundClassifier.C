#include <iostream>
using namespace std;
#include <vector>

#include "FinalStatesClass.C"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

// classifying singals, and storing the truth level indices
Int_t ClassifySignal(TClonesArray* branchParticle, TLorentzVector* BTrue, TLorentzVector* HcTrue, TLorentzVector* muTrue, TVector3* v3HcTrue, TVector3* v3MuTrue, TVector3* v3BTrue, iFinalStates* iFSTrue, Int_t mode) {
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
                if (particle0->M1 == -1) continue;
                particle0M = (GenParticle*)branchParticle->At(particle0->M1);
                if (particle0M->M1 == -1) continue;
                particle0MM = (GenParticle*)branchParticle->At(particle0M->M1);

                if (abs(particle0->PID) == 431 && particle0->M1 == iBTrue) {
                    iCTrue_i = ip2;
                    foundCTrue_i = 1;
                }
                if (abs(particle0->PID) == 333 && abs(particle0M->PID) == 431 && particle0M->M1 == iBTrue) {
                    iC0True_i = ip2;
                    foundC0True_i = 1;
                }
                if (particle0->PID == -321 && abs(particle0M->PID) == 333 && abs(particle0MM->PID) == 431 && particle0MM->M1 == iBTrue) {
                    iC1True_i = ip2;
                    foundC1True_i = 1;
                }
                if (particle0->PID == 321 && abs(particle0M->PID) == 333 && abs(particle0MM->PID) == 431 && particle0MM->M1 == iBTrue) {
                    iC2True_i = ip2;
                    foundC2True_i = 1;
                }
                if (abs(particle0->PID) == 211 && abs(particle0M->PID) == 431 && particle0M->M1 == iBTrue) {
                    iC3True_i = ip2;
                    foundC3True_i = 1;
                }

                if (mode == 1 && (abs(particle0->PID) == 13 && abs(particle0M->PID) == 15 && particle0M->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
                if (mode == 2 && (abs(particle0->PID) == 13 && particle0->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
            }
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
        TVector3 v3HcTrue_, v3MuTrue_, v3BTrue_;
        v3HcTrue_.SetXYZ(particleC1->X, particleC1->Y, particleC1->Z);
        v3MuTrue_.SetXYZ(particleMu->X, particleMu->Y, particleMu->Z);
        v3BTrue_.SetXYZ(particleHc->X, particleHc->Y, particleHc->Z);
        // cout << " X: " << particleHc->X << "\n";

        *iFSTrue = iFSTrue_;
        *BTrue = BTrue_;
        *HcTrue = HcTrue_;
        *muTrue = muTrue_;
        *v3HcTrue = v3HcTrue_;
        *v3MuTrue = v3MuTrue_;
        *v3BTrue = v3BTrue_;
        return 1;
    } else {
        return 0;
    }
}

// classifying excited singals (D_s^*), and storing the truth level indices
Int_t ClassifyExcitedSignal(TClonesArray* branchParticle, TLorentzVector* BTrue, TLorentzVector* HcTrue, TLorentzVector* muTrue, TLorentzVector* phoTrue, TVector3* v3HcTrue, TVector3* v3MuTrue, TVector3* v3BTrue, iFinalStates* iFSTrue, Int_t mode) {
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
                if (particle0->M1 == -1) continue;
                particle0M = (GenParticle*)branchParticle->At(particle0->M1);
                if (particle0M->M1 == -1) continue;
                particle0MM = (GenParticle*)branchParticle->At(particle0M->M1);
                if (particle0MM->M1 == -1) continue;
                particle0MMM = (GenParticle*)branchParticle->At(particle0MM->M1);

                if (abs(particle0->PID) == 431 && abs(particle0M->PID) == 433 && particle0M->M1 == iBTrue) {
                    iCTrue_i = ip2;
                    foundCTrue_i = 1;
                }
                if (abs(particle0->PID) == 333 && abs(particle0->PID) == 431 && abs(particle0M->PID) == 433 && particle0MM->M1 == iBTrue) {
                    iC0True_i = ip2;
                    foundC0True_i = 1;
                }
                if (particle0->PID == -321 && abs(particle0M->PID) == 333 && abs(particle0MM->PID) == 431 && abs(particle0MMM->PID) == 433 && particle0MMM->M1 == iBTrue) {
                    iC1True_i = ip2;
                    foundC1True_i = 1;
                }
                if (particle0->PID == 321 && abs(particle0M->PID) == 333 && abs(particle0MM->PID) == 431 && abs(particle0MMM->PID) == 433 && particle0MMM->M1 == iBTrue) {
                    iC2True_i = ip2;
                    foundC2True_i = 1;
                }
                if (abs(particle0->PID) == 211 && abs(particle0M->PID) == 431 && abs(particle0MM->PID) == 433 && particle0MM->M1 == iBTrue) {
                    iC3True_i = ip2;
                    foundC3True_i = 1;
                }
                if (mode == 1 && (abs(particle0->PID) == 13 && abs(particle0M->PID) == 15 && particle0M->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
                if (mode == 2 && (abs(particle0->PID) == 13 && particle0->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
            }
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
        if (abs(particlePho->PID) != 22) continue;
        if (particlePho->M1 == particleC->M1) {
            foundPhoTrue = 1;
            iPhoTrue = ipho;
        }
    }

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
        TVector3 v3HcTrue_, v3MuTrue_, v3BTrue_;
        v3HcTrue_.SetXYZ(particleC1->X, particleC1->Y, particleC1->Z);
        v3MuTrue_.SetXYZ(particleMu->X, particleMu->Y, particleMu->Z);
        v3BTrue_.SetXYZ(particlePho->X, particlePho->Y, particlePho->Z);
        // cout << " particlePho->X: " << particlePho->X << "\n";
        // cout << " particleHc->X: " << particleHc->X << "\n";
        // cout << endl;

        *iFSTrue = iFSTrue_;
        *BTrue = BTrue_;
        *HcTrue = HcTrue_;
        *muTrue = muTrue_;
        *phoTrue = phoTrue_;
        *v3HcTrue = v3HcTrue_;
        *v3MuTrue = v3MuTrue_;
        *v3BTrue = v3BTrue_;
        return 1;
    } else {
        return 0;
    }
}

// classify if the event is the targeted bkg (from inclusive samples)
Int_t ClassifyBkg(TClonesArray* branchParticle, const string type, Int_t* bkgHaveFromDsstar, vector<TLorentzVector>* bkgPhotons) {
    // Inclusive samepls (for Comb+Cascade, Comb, Casde, Inclusive)
    GenParticle* particle;
    Int_t nParticles = branchParticle->GetEntries();

    //=== check c-hadron ===
    Int_t nPosHc = 0;  // number of positive targeted c-hadron in the truth level (Jpsi, Ds, Lambdac)
    Int_t nNegHc = 0;  // number of negative targeted c-hadron in the truth level (Jpsi, Ds, Lambdac)
    for (int ip = 0; ip < nParticles; ip++) {
        particle = (GenParticle*)branchParticle->At(ip);
        if (particle->PID == 431) nPosHc += 1;   // check number of pos lambdac
        if (particle->PID == -431) nNegHc += 1;  // check number of neg lambdac
    }
    if (not(nPosHc + nNegHc >= 1)) return 0;  // must have at least one lambdac

    //=== check muons ===
    Int_t nPosMu = 0;  // number of positive muons
    Int_t nNegMu = 0;  // number of negative muons
    for (int ip = 0; ip < nParticles; ip++) {
        particle = (GenParticle*)branchParticle->At(ip);
        if (particle->PID == -13) nPosMu += 1;
        if (particle->PID == 13) nNegMu += 1;
    }
    // total must be at least 3, (2 from Jpsi, 1 unapired), and at least 1 pos,
    // 1 neg from Jpsi
    if (not((nPosHc > 0 && nNegMu > 0) || (nNegHc > 0 && nPosMu > 0))) return 0;

    //=== check the decay chain ===
    Int_t isLeptonFromB = 0;      // check if the lepton is from b-hadron
    Int_t isLeptonSemiFromB = 0;  // check if the lepton is semileptonically decay from b-hadron
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
            if (abs(particle1M->PID) == 15) iLepton = particle1->M1;
            particleL = (GenParticle*)branchParticle->At(iLepton);
            Int_t iLeptonMother = particleL->M1;
            particle1M = (GenParticle*)branchParticle->At(iLeptonMother);

            // finding the b-hadron that lepton from
            Int_t BHadron_idx = 99999;
            nLeptonStepFromB = 0;
            while (true) {
                if (int(abs(particle1M->PID) / 100) == 5 ||
                    int(abs(particle1M->PID) / 1000) == 5 ||
                    int(abs(particle1M->PID) / 10000) == 5 ||
                    int((abs(particle1M->PID) % 1000) / 100) == 5) {
                    BHadron_idx = iLeptonMother;
                    break;
                } else if (int(particle1M->PID / 10) == 0) {
                    break;
                }
                nLeptonStepFromB += 1;
                iLeptonMother = particle1M->M1;
                particle1M = (GenParticle*)branchParticle->At(iLeptonMother);
            }
            if (BHadron_idx != 99999) isLeptonFromB = 1;
            if (BHadron_idx != 99999 && nLeptonStepFromB == 0) isLeptonSemiFromB = 1;

            // if lepton is from b-hadron, then find if c-hadron is from the
            // same mother
            if (isLeptonFromB == 1) {
                for (int ip2 = 0; ip2 < nParticles; ip2++) {
                    particle2 = (GenParticle*)branchParticle->At(ip2);
                    if (abs(particle2->PID) == 431) {
                        Int_t iC = ip2;
                        particleC = (GenParticle*)branchParticle->At(iC);
                        Int_t iCMother = particleC->M1;
                        particle2M = (GenParticle*)branchParticle->At(iCMother);
                        Int_t iCMotherMother = particle2M->M1;
                        if (particle2M->M1 == -1) continue;
                        // find if the bkg has Ds* or just Ds
                        // cout << " particleCM: " << particle2M->PID << "\n";
                        if (abs(particle2M->PID) == 433) {
                            *bkgHaveFromDsstar = 1;
                        } else {
                            *bkgHaveFromDsstar = 0;
                        }

                        particle2MM = (GenParticle*)branchParticle->At(particle2M->M1);

                        if (iCMother == BHadron_idx && isLeptonSemiFromB == 1) {
                            // testing
                            // Int_t nParticlesFromB = 0;
                            // GenParticle* testingBHadron = (GenParticle*)branchParticle->At(BHadron_idx);
                            // if (abs(testingBHadron->PID) != 531) cout << " text: " << testingBHadron->PID << "\n";
                            // for (Int_t ifromB = 0; ifromB < nParticles; ifromB++) {
                            // GenParticle* particleFromB = (GenParticle*)branchParticle->At(ifromB);
                            // if (particleFromB->M1 == BHadron_idx) nParticlesFromB += 1;
                            //}
                            // if (nParticlesFromB != 3) cout << " nParticlesFromB: " << nParticlesFromB << "\n";
                            // end testing
                            isSignal = 1;
                            continue;
                        }
                        if (abs(particle2M->PID) == 433 &&
                            iCMotherMother == BHadron_idx &&
                            isLeptonSemiFromB == 1) {
                            isSignal = 1;
                            continue;
                        }

                        while (true) {
                            if (iCMother == BHadron_idx) {
                                isCFromB = 1;
                                break;
                            }
                            if (int(particle2M->PID / 10) == 0) {
                                break;
                            }
                            iCMother = particle2M->M1;
                            particle2M = (GenParticle*)branchParticle->At(iCMother);
                        }
                        if (isCFromB == 1 || isSignal == 1) break;
                    }
                }
                if (isCFromB == 1 || isSignal == 1) break;
            }
        }
    }
    if (isSignal == 1) return 0;

    vector<TLorentzVector> storePhotons;
    // cout << " *bkgHaveFromDsstar: " << *bkgHaveFromDsstar << "\n";
    if (*bkgHaveFromDsstar == 1) {
        // cout << " *bkgHaveFromDsstar: " << *bkgHaveFromDsstar << "\n";
        Int_t nParticles = branchParticle->GetEntries();
        for (Int_t iph = 0; iph < nParticles; iph++) {
            GenParticle* photon = (GenParticle*)branchParticle->At(iph);
            if (abs(photon->PID) != 22) continue;
            GenParticle* photonM = (GenParticle*)branchParticle->At(photon->M1);
            // cout << " photonM->PID: " << photonM->PID << "\n";
            if (abs(photonM->PID) == 433) {
                TLorentzVector pho;
                pho.SetPtEtaPhiE(photon->PT, photon->Eta, photon->Phi, photon->E);
                // cout << " found : \n";
                storePhotons.push_back(pho);
            }
            // cout << " photon->M1: " << photon->M1 << "\n";
        }
        if (storePhotons.size() != 0) {
            *bkgPhotons = storePhotons;
            // cout << " storePhotn0: " << storePhotons[0].E() << "\n";
        }
    }
    // cout << " c pid: " << particleC->PID << "\n";
    // cout << " m1 index: " << particleC->M1 << "\n";
    // GenParticle* particleCM = (GenParticle*)branchParticle->At(particleC->M1);
    // cout << " c m pid: " << particleCM->PID << "\n";
    // if (abs(particleCM->PID) == 433) {
    //*bkgHaveFromDsstar = 1;
    //} else {
    //*bkgHaveFromDsstar = 0;
    //}
    Int_t isComb = 0;
    Int_t isCascade = 0;
    Int_t isInclusive = 0;

    if (isLeptonSemiFromB == 1 && isCFromB == 1) isInclusive = 1;
    if (isLeptonFromB == 1 && isLeptonSemiFromB == 0 && isCFromB == 1) isCascade = 1;
    if (isCFromB == 0) isComb = 1;
    if (isInclusive + isCascade + isComb != 1) cout << "HAVE BUG IN CLASSIFYING BKG!";

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

Int_t ClassifyMisID(TClonesArray* branchParticle, Int_t* nPi_, Int_t* bkgHaveFromDsstar) {
    // skip signal events
    TLorentzVector dummy;
    TVector3 dummy2;
    iFinalStates dummy3;
    if (ClassifySignal(branchParticle, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy2, &dummy3, 1) == 1) {
        return 0;
    } else if (ClassifySignal(branchParticle, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy2, &dummy3, 2) == 1) {
        return 0;
    } else if (ClassifyExcitedSignal(branchParticle, &dummy, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy2, &dummy3, 1) == 1) {
        return 0;
    } else if (ClassifyExcitedSignal(branchParticle, &dummy, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy2, &dummy3, 2) == 1) {
        return 0;
    } else {
        GenParticle* particle;
        Int_t nParticles = branchParticle->GetEntries();

        //=== check c-hadron ===
        Int_t nHc = 0;  // number of targeted c-hadron in the truth level (Jpsi, Ds, Lambdac)
        *bkgHaveFromDsstar = 0;
        for (int ip = 0; ip < nParticles; ip++) {
            particle = (GenParticle*)branchParticle->At(ip);
            if (abs(particle->PID) == 431) {
                nHc += 1;  // check number of Jpsi
                GenParticle* particle2M = (GenParticle*)branchParticle->At(particle->M1);
                if (abs(particle2M->PID) == 433) *bkgHaveFromDsstar = 1;
            }
        }
        if (not(nHc >= 1)) return 0;  // must have at least one Jpsi

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
        if (not(nPi >= 2)) return 0;

        *nPi_ = nPi;
        return 1;
    }
}

