
#include <iostream>
using namespace std;
#include "FinalStatesClass.C"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

Int_t ClassifySignal(TClonesArray *branchParticle, TLorentzVector *BTrue, TLorentzVector *HcTrue, TLorentzVector *muTrue, TVector3 *v3HcTrue, TVector3 *v3MuTrue, iFinalStates *iFSTrue, Int_t mode) {
    GenParticle *particle;
    GenParticle *particleB;
    GenParticle *particleC;
    GenParticle *particle0;
    GenParticle *particle0M;
    GenParticle *particle0MM;
    Int_t nParticles = branchParticle->GetEntries();

    Int_t iBTrue, foundBTrue = 0;
    Int_t iCTrue, foundCTrue = 0;
    Int_t iC1True, foundC1True = 0;
    Int_t iC2True, foundC2True = 0;
    Int_t iC3True, foundC3True = 0;
    Int_t iMuTrue, foundMuTrue = 0;
    for (int ipB = 0; ipB < nParticles; ipB++) {
        foundBTrue = 0;
        particleB = (GenParticle *)branchParticle->At(ipB);
        if (abs(particleB->PID) == 5122) {
            iBTrue = ipB;
            foundBTrue = 1;

            Int_t iCTrue_i, foundCTrue_i = 0;
            Int_t iC1True_i, foundC1True_i = 0;
            Int_t iC2True_i, foundC2True_i = 0;
            Int_t iC3True_i, foundC3True_i = 0;
            Int_t iMuTrue_i, foundMuTrue_i = 0;
            for (int ip2 = 0; ip2 < nParticles; ip2++) {
                particle0 = (GenParticle *)branchParticle->At(ip2);
                if (particle0->M1 == -1) {
                    continue;
                }
                particle0M = (GenParticle *)branchParticle->At(particle0->M1);

                if (abs(particle0->PID) == 4122 && particle0->M1 == iBTrue) {
                    iCTrue_i = ip2;
                    foundCTrue_i = 1;
                }
                if (abs(particle0->PID) == 2212 && abs(particle0M->PID) == 4122 && particle0M->M1 == iBTrue) {
                    iC1True_i = ip2;
                    foundC1True_i = 1;
                }
                if (abs(particle0->PID) == 321 && abs(particle0M->PID) == 4122 && particle0M->M1 == iBTrue) {
                    iC2True_i = ip2;
                    foundC2True_i = 1;
                }
                if (abs(particle0->PID) == 211 && abs(particle0M->PID) == 4122 && particle0M->M1 == iBTrue) {
                    iC3True_i = ip2;
                    foundC3True_i = 1;
                }

                if (particle0M->M1 == -1) {
                    continue;
                }
                particle0MM = (GenParticle *)branchParticle->At(particle0M->M1);
                if (abs(particle0->PID) == 2212 && abs(particle0MM->PID) == 4122 && particle0MM->M1 == iBTrue) {
                    iC1True_i = ip2;
                    foundC1True_i = 1;
                }
                if (abs(particle0->PID) == 321 && abs(particle0MM->PID) == 4122 && particle0MM->M1 == iBTrue) {
                    iC2True_i = ip2;
                    foundC2True_i = 1;
                }
                if (abs(particle0->PID) == 211 && abs(particle0MM->PID) == 4122 && particle0MM->M1 == iBTrue) {
                    iC3True_i = ip2;
                    foundC3True_i = 1;
                }

                /*
                if ((abs(particle0->PID) == 13 && particle0->M1 == iBTrue) ||
                    (abs(particle0->PID) == 13 && abs(particle0M->PID) == 15 && particle0M->M1 == iBTrue)) {
                    iMuTrue_i = ip2; foundMuTrue_i = 1;
                }
              */

                if (mode == 1 && (abs(particle0->PID) == 13 && abs(particle0M->PID) == 15 && particle0M->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
                if (mode == 2 && (abs(particle0->PID) == 13 && particle0->M1 == iBTrue)) {
                    iMuTrue_i = ip2;
                    foundMuTrue_i = 1;
                }
            }
            // cout <<  foundCTrue_i <<"; "<< foundC1True_i <<"; "<< foundC2True_i <<"; "<< foundC3True_i<< "; " << foundMuTrue_i<<endl;
            if (foundCTrue_i == 1 && foundC1True_i == 1 && foundC2True_i == 1 && foundC3True_i == 1 && foundMuTrue_i == 1) {
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
    // cout <<  foundBTrue <<"; "<< foundCTrue <<"; "<< foundCTrue <<"; "<< foundC2True <<"; "<< foundC3True << "; " << foundMuTrue <<endl;
    if (foundBTrue == 1 && foundCTrue == 1 && foundCTrue == 1 && foundC2True == 1 && foundC3True == 1 && foundMuTrue == 1) {
        iFinalStates iFSTrue_;
        iFSTrue_.iP = iC1True;
        iFSTrue_.iK = iC2True;
        iFSTrue_.iPi = iC3True;
        iFSTrue_.iMu = iMuTrue;

        TLorentzVector BTrue_, HcTrue_, muTrue_;
        GenParticle *particleHc;
        GenParticle *particleMu;
        GenParticle *particleC1;
        particleB = (GenParticle *)branchParticle->At(iBTrue);
        BTrue_.SetPtEtaPhiE(particleB->PT, particleB->Eta, particleB->Phi, particleB->E);
        particleHc = (GenParticle *)branchParticle->At(iCTrue);
        HcTrue_.SetPtEtaPhiE(particleHc->PT, particleHc->Eta, particleHc->Phi, particleHc->E);
        particleMu = (GenParticle *)branchParticle->At(iMuTrue);
        muTrue_.SetPtEtaPhiE(particleMu->PT, particleMu->Eta, particleMu->Phi, particleMu->E);

        particleC1 = (GenParticle *)branchParticle->At(iC1True);
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

//---------- return if the event is the targeted bkg (from inclusive samples) {{{
Int_t ClassifyBkg(TClonesArray *branchParticle, const string type) {
    // Inclusive samepls (for Comb+Cascade, Comb, Casde, Inclusive)
    GenParticle *particle;
    Int_t nParticles = branchParticle->GetEntries();

    //=== check c-hadron ===
    Int_t nPosHc = 0;  // number of positive targeted c-hadron in the truth level (Jpsi, Ds, Lambdac)
    Int_t nNegHc = 0;  // number of negative targeted c-hadron in the truth level (Jpsi, Ds, Lambdac)
    for (int ip = 0; ip < nParticles; ip++) {
        particle = (GenParticle *)branchParticle->At(ip);
        if (particle->PID == 4122) {
            nPosHc += 1;
        }  // check number of pos lambdac
        if (particle->PID == -4122) {
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
        particle = (GenParticle *)branchParticle->At(ip);
        if (particle->PID == -13) {
            nPosMu += 1;
        }
        if (particle->PID == 13) {
            nNegMu += 1;
        }
    }
    // total must be at least 3, (2 from Jpsi, 1 unapired), and at least 1 pos, 1 neg from Jpsi
    if (not((nPosHc > 0 && nNegMu > 0) || (nNegHc > 0 && nPosMu > 0))) {
        return 0;
    }

    //=== check the decay chain ===
    Int_t isLeptonFromB = 0;      // check if the lepton is from b-hadron
    Int_t isLeptonSemiFromB = 0;  // check if the lepton is semileptonically decay from b-hadron
    Int_t nLeptonStepFromB = 0;
    Int_t isCFromB = 0;  // check if c-hadron is from B
    Int_t isSignal = 0;  // chekc if event has signal in it

    GenParticle *particle1;
    GenParticle *particle1M;  // particle1's mother
    GenParticle *particle2;
    GenParticle *particle2M;  // particle2's mother
    GenParticle *particleL;   // lepton
    GenParticle *particleC;   // c-hadron
    for (int ip1 = 0; ip1 < nParticles; ip1++) {
        particle1 = (GenParticle *)branchParticle->At(ip1);
        // identify the muon
        if (abs(particle1->PID) == 13) {
            isLeptonFromB = 0;
            isLeptonSemiFromB = 0;
            isCFromB = 0;

            // mother of particle0
            particle1M = (GenParticle *)branchParticle->At(particle1->M1);

            // identify the lepton (mu/tau)
            Int_t iLepton = ip1;
            if (abs(particle1M->PID) == 15) {
                iLepton = particle1->M1;
            }
            particleL = (GenParticle *)branchParticle->At(iLepton);
            Int_t iLeptonMother = particleL->M1;
            particle1M = (GenParticle *)branchParticle->At(iLeptonMother);
            // cout << "Lepton: " << particleL->PID << " (" << iLepton << ") from: ";

            // finding the b-hadron that lepton from
            Int_t BHadron_idx = 99999;
            nLeptonStepFromB = 0;
            while (true) {
                // cout << particle1M->PID << " < ";
                if (int(abs(particle1M->PID) / 100) == 5 || int(abs(particle1M->PID) / 1000) == 5 || int(abs(particle1M->PID) / 10000) == 5 ||
                    int((abs(particle1M->PID) % 1000) / 100) == 5) {
                    BHadron_idx = iLeptonMother;
                    // cout << "\n From b-hadron (IDX): " << particle1M->PID << " (" << BHadron_idx << ")" << endl;
                    break;
                } else if (int(particle1M->PID / 10) == 0) {
                    // cout << "\nNo b-hadron, terminated at: " << particle1M->PID << endl;
                    break;
                }
                nLeptonStepFromB += 1;
                iLeptonMother = particle1M->M1;
                particle1M = (GenParticle *)branchParticle->At(iLeptonMother);
            }
            if (BHadron_idx != 99999) {
                isLeptonFromB = 1;
            }
            if (BHadron_idx != 99999 && nLeptonStepFromB == 0) {
                isLeptonSemiFromB = 1;
            }

            // if lepton is from b-hadron, then find if c-hadron is from the same mother
            if (isLeptonFromB == 1) {
                for (int ip2 = 0; ip2 < nParticles; ip2++) {
                    particle2 = (GenParticle *)branchParticle->At(ip2);
                    if (abs(particle2->PID) == 4122) {
                        Int_t iC = ip2;
                        particleC = (GenParticle *)branchParticle->At(iC);
                        Int_t iCMother = particleC->M1;
                        particle2M = (GenParticle *)branchParticle->At(iCMother);

                        if (iCMother == BHadron_idx && isLeptonSemiFromB == 1) {
                            isSignal = 1;
                            continue;
                        }
                        // cout << "c-hadron " << particleC->PID << " (" << iCMother << ") from: ";

                        while (true) {
                            // cout << particle2M->PID << " (" << iCMother << ")" << " < ";
                            if (iCMother == BHadron_idx) {
                                // cout << "/nc-hadron from the smae b-hadron: " << particle2M->PID << endl;
                                isCFromB = 1;
                                break;
                            }
                            if (int(particle2M->PID / 10) == 0) {
                                // cout << "\nc-hadron not from same b-hadron terminate at: " << particle2M->PID << endl;
                                break;
                            }
                            iCMother = particle2M->M1;
                            particle2M = (GenParticle *)branchParticle->At(iCMother);
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

    // cout << "Signal: " << isSignal << "; Inclsuive: " << isInclusive << "; Cascade: " << isCascade << "; Comb:" << isComb << endl;
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

//---------- return if the event is misID bkg {{{
Int_t ClassifyMisID(TClonesArray *branchParticle) {
    // skip signal events
    TLorentzVector dummy;
    TVector3 dummy2;
    iFinalStates dummy3;
    if (ClassifySignal(branchParticle, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy3, 1) == 1) {
        return 0;
    } else if (ClassifySignal(branchParticle, &dummy, &dummy, &dummy, &dummy2, &dummy2, &dummy3, 2) == 1) {
        return 0;
    } else {
        GenParticle *particle;
        Int_t nParticles = branchParticle->GetEntries();

        //=== check c-hadron ===
        Int_t nHc = 0;  // number of targeted c-hadron in the truth level (Jpsi, Ds, Lambdac)
        for (int ip = 0; ip < nParticles; ip++) {
            particle = (GenParticle *)branchParticle->At(ip);
            if (abs(particle->PID) == 4122) {
                nHc += 1;
            }  // check number of Jpsi
        }
        if (not(nHc >= 1)) {
            return 0;
        }  // must have at least one Jpsi

        //=== check pions ===
        Int_t nPi = 0;  // number of pions (displaced)
        for (int ip = 0; ip < nParticles; ip++) {
            particle = (GenParticle *)branchParticle->At(ip);
            if (abs(particle->PID) == 211 && (not(particle->X == 0 && particle->Y == 0 && particle->Z == 0))) {
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

