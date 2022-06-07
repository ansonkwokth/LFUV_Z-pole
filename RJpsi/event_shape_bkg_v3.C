#include <cmath>
#include <iostream>

using namespace std;

void event_shape_bkg_v3() {
    // TFile* f = new TFile("./Bcjpsitaunu_50m.root");
    // TFile* f = new TFile("./BcJpsimunu0-2.root");
    TFile* f = new TFile("./RJpsi_comb_200m_seed1.root");
    TTree* tree = (TTree*)f->Get("Delphes");

    // TFile fea("./features/2b_FWmoments.root", "recreate");
    // TFile fea("./features/2b2c_FWmoments.root", "recreate");
    // TFile fea("./features/4b_FWmoments.root", "recreate");
    TTree tr("t", "features");
    Float_t H_EE0_store, H_EE1_store, H_EE2_store, H_EE3_store, H_EE4_store, H_EE5_store;
    Float_t H_EE6_store, H_EE7_store, H_EE8_store, H_EE9_store, H_EE10_store;
    tr.Branch("H_EE0", &H_EE0_store, "H_EE0/F");
    tr.Branch("H_EE1", &H_EE1_store, "H_EE1/F");
    tr.Branch("H_EE2", &H_EE2_store, "H_EE2/F");
    tr.Branch("H_EE3", &H_EE3_store, "H_EE3/F");
    tr.Branch("H_EE4", &H_EE4_store, "H_EE4/F");
    tr.Branch("H_EE5", &H_EE5_store, "H_EE5/F");
    tr.Branch("H_EE6", &H_EE6_store, "H_EE6/F");
    tr.Branch("H_EE7", &H_EE7_store, "H_EE7/F");
    tr.Branch("H_EE8", &H_EE8_store, "H_EE8/F");
    tr.Branch("H_EE9", &H_EE9_store, "H_EE9/F");
    tr.Branch("H_EE10", &H_EE10_store, "H_EE10/F");

    // TFile fea2("./features/angles.root", "recreate");
    // TFile fea2("./features/angles_v2.root", "recreate");
    TTree tr2("t", "features");
    Int_t have_Bc_store;
    Int_t num_b_store, num_c_store, env_idx_store;
    Float_t angle1_store, angle2_store, E1_store, E2_store;
    tr2.Branch("env_idx", &env_idx_store, "env_idx/I");
    tr2.Branch("have_Bc", &have_Bc_store, "have_Bc/I");
    tr2.Branch("num_b", &num_b_store, "num_b/I");
    tr2.Branch("num_c", &num_c_store, "num_c/I");
    tr2.Branch("angle1", &angle1_store, "angle1/F");
    tr2.Branch("angle2", &angle2_store, "angle2/F");
    tr2.Branch("E1", &E1_store, "E1/F");
    tr2.Branch("E2", &E2_store, "E2/F");
    tr2.Branch("H_EE0", &H_EE0_store, "H_EE0/F");
    tr2.Branch("H_EE1", &H_EE1_store, "H_EE1/F");
    tr2.Branch("H_EE2", &H_EE2_store, "H_EE2/F");
    tr2.Branch("H_EE3", &H_EE3_store, "H_EE3/F");
    tr2.Branch("H_EE4", &H_EE4_store, "H_EE4/F");
    tr2.Branch("H_EE5", &H_EE5_store, "H_EE5/F");
    tr2.Branch("H_EE6", &H_EE6_store, "H_EE6/F");
    tr2.Branch("H_EE7", &H_EE7_store, "H_EE7/F");
    tr2.Branch("H_EE8", &H_EE8_store, "H_EE8/F");
    tr2.Branch("H_EE9", &H_EE9_store, "H_EE9/F");
    tr2.Branch("H_EE10", &H_EE10_store, "H_EE10/F");

    Int_t nEntries;
    nEntries = tree->GetEntries();

    // load information
    const int max_val = 1000;
    Float_t m_K0 = 0.497661;

    // track information
    Int_t tr_size, tr_pid[max_val], tr_charge[max_val];
    Float_t tr_X[max_val], tr_Y[max_val], tr_Z[max_val], tr_T[max_val], tr_pt[max_val], tr_eta[max_val], tr_phi[max_val], tr_p[max_val];

    tree->SetBranchAddress("Track_size", &tr_size);
    tree->SetBranchAddress("Track.PID", tr_pid);
    tree->SetBranchAddress("Track.X", tr_X);
    tree->SetBranchAddress("Track.Y", tr_Y);
    tree->SetBranchAddress("Track.Z", tr_Z);
    tree->SetBranchAddress("Track.T", tr_T);
    tree->SetBranchAddress("Track.Charge", tr_charge);
    tree->SetBranchAddress("Track.PT", tr_pt);
    tree->SetBranchAddress("Track.Eta", tr_eta);
    tree->SetBranchAddress("Track.Phi", tr_phi);
    tree->SetBranchAddress("Track.P", tr_p);

    // photon
    Int_t efp_size;
    Float_t efp_et[max_val], efp_eta[max_val], efp_phi[max_val], efp_E[max_val];
    tree->SetBranchAddress("EFlowPhoton_size", &efp_size);
    tree->SetBranchAddress("EFlowPhoton.ET", efp_et);
    tree->SetBranchAddress("EFlowPhoton.Eta", efp_eta);
    tree->SetBranchAddress("EFlowPhoton.Phi", efp_phi);
    tree->SetBranchAddress("EFlowPhoton.E", efp_E);

    // neutral hadron
    Int_t efnh_size;
    Float_t efnh_et[max_val], efnh_eta[max_val], efnh_phi[max_val], efnh_E[max_val];
    tree->SetBranchAddress("EFlowNeutralHadron_size", &efnh_size);
    tree->SetBranchAddress("EFlowNeutralHadron.ET", efnh_et);
    tree->SetBranchAddress("EFlowNeutralHadron.Eta", efnh_eta);
    tree->SetBranchAddress("EFlowNeutralHadron.Phi", efnh_phi);
    tree->SetBranchAddress("EFlowNeutralHadron.E", efnh_E);

    // particle
    Int_t p_size, p_pid[max_val], p_M1[max_val];
    Float_t p_X[max_val], p_Y[max_val], p_Z[max_val], p_pt[max_val], p_eta[max_val], p_phi[max_val], p_E[max_val];
    tree->SetBranchAddress("Particle_size", &p_size);
    tree->SetBranchAddress("Particle.PID", p_pid);
    tree->SetBranchAddress("Particle.X", p_X);
    tree->SetBranchAddress("Particle.Y", p_Y);
    tree->SetBranchAddress("Particle.Z", p_Z);
    tree->SetBranchAddress("Particle.PT", p_pt);
    tree->SetBranchAddress("Particle.Eta", p_eta);
    tree->SetBranchAddress("Particle.Phi", p_phi);
    tree->SetBranchAddress("Particle.E", p_E);
    tree->SetBranchAddress("Particle.M1", p_M1);

    // TCanvas* c = new TCanvas();
    // c->Divide(2, 2, 0.0000001, 0.0000001);
    // TH1F* hist_H_EE_1 = new TH1F("H_EE_1", "H_{EE,1}", 24, -0.1, 1.1);
    // TH1F* hist_H_EE_2 = new TH1F("H_EE_2", "H_{EE,2}", 24, -0.1, 1.1);
    // TH1F* hist_H_EE_3 = new TH1F("H_EE_3", "H_{EE,3}", 24, -0.1, 1.1);
    // TH1F* hist_H_EE_4 = new TH1F("H_EE_4", "H_{EE,4}", 24, -0.1, 1.1);

    int event_count = 0;
    int look_back = 7;

    //##########################################
    //##########    looping events    ##########
    //##########################################
    // for (int i_en = 0; i_en < nEntries; i_en++) {
    for (int i_en = 0; i_en < 100; i_en++) {
        tree->GetEntry(i_en);

        if (i_en % 10000 == 0) {
            cout << i_en << "/" << nEntries << endl;
        }

        // exclude signals
        int truth_Bc = 0;
        int truth_from_Jpsi = 0;
        int truth_from_unpaired = 0;
        int have_mu = 0;
        for (int ii = 0; ii < p_size; ii++) {
            if (abs(p_pid[ii]) == 13 && abs(p_pid[p_M1[ii]]) == 443 && abs(p_pid[p_M1[p_M1[ii]]]) == 541) {
                truth_from_Jpsi += 1;
            }
            // if (abs(p_pid[ii]) == 13 && abs(p_pid[p_M1[ii]]) == 443 && (abs(p_pid[p_M1[p_M1[ii]]]) == 10441 || abs(p_pid[p_M1[p_M1[ii]]]) == 20443 || abs(p_pid[p_M1[p_M1[ii]]]) == 445) && abs(p_pid[p_M1[p_M1[p_M1[ii]]]]) == 541) { truth_from_Jpsi += 1; }
            if ((abs(p_pid[ii]) == 13 && abs(p_pid[p_M1[ii]]) == 541) ||
                (abs(p_pid[ii]) == 13 && abs(p_pid[p_M1[ii]]) == 15 && abs(p_pid[p_M1[p_M1[ii]]]) == 541)) {
                truth_from_unpaired += 1;
            }
            if (abs(p_pid[ii]) == 541) {
                truth_Bc = 1;
            }
            if (abs(p_pid[ii]) == 13) {
                have_mu = 1;
            }
        }
        if (truth_from_Jpsi == 2 && truth_from_unpaired == 1) {
            continue;
        }
        have_Bc_store = truth_Bc;

        // exclude excited states
        truth_Bc = 0;
        truth_from_Jpsi = 0;
        truth_from_unpaired = 0;
        have_mu = 0;
        for (int ii = 0; ii < p_size; ii++) {
            // if (abs(p_pid[ii]) == 13 && abs(p_pid[p_M1[ii]]) == 443 && abs(p_pid[p_M1[p_M1[ii]]]) == 541) { truth_from_Jpsi += 1; }
            if (abs(p_pid[ii]) == 13 && abs(p_pid[p_M1[ii]]) == 443 && (abs(p_pid[p_M1[p_M1[ii]]]) == 10441 || abs(p_pid[p_M1[p_M1[ii]]]) == 20443 || abs(p_pid[p_M1[p_M1[ii]]]) == 445) && abs(p_pid[p_M1[p_M1[p_M1[ii]]]]) == 541) {
                truth_from_Jpsi += 1;
            }
            if ((abs(p_pid[ii]) == 13 && abs(p_pid[p_M1[ii]]) == 541) ||
                (abs(p_pid[ii]) == 13 && abs(p_pid[p_M1[ii]]) == 15 && abs(p_pid[p_M1[p_M1[ii]]]) == 541)) {
                truth_from_unpaired += 1;
            }
            if (abs(p_pid[ii]) == 541) {
                truth_Bc = 1;
            }
            if (abs(p_pid[ii]) == 13) {
                have_mu = 1;
            }
        }
        if (truth_from_Jpsi == 2 && truth_from_unpaired == 1) {
            continue;
        }

        // cout << truth_from_Jpsi << ": " << truth_from_unpaired << endl;
        // if (truth_from_Jpsi != 2) { continue; }
        // if (truth_from_unpaired != 1) { continue; }
        // cout <<"............" << have_mu << endl;
        have_Bc_store = truth_Bc;

        /*
        int c_size = 0;
        int jpsi_found = 0, B0_size = 0, Bpm_size = 0, Bs0_size = 0, Bc_size = 0, Lambda_size = 0, Xi_size = 0;

        int c1_idx = 0, c2_idx = 0;

        int background_found = 0;
        for (int p_i = 0; p_i < p_size; p_i++) {
                // ***************************** counting B ***********************************
                if (abs(p_pid[p_i]) == 511) { B0_size += 1; }
                if (abs(p_pid[p_i]) == 521) { Bpm_size += 1; }
                if (abs(p_pid[p_i]) == 531) { Bs0_size += 1; }
                if (abs(p_pid[p_i]) == 541) { Bc_size += 1; }
                if (abs(p_pid[p_i]) == 5122) { Lambda_size += 1; }
                if (abs(p_pid[p_i]) == 5132) { Xi_size += 1; }
                if (abs(p_pid[p_i]) == 5232) { Xi_size += 1; }
                // ******************************** counting c *************************************
                if (abs(p_pid[p_M1[p_i]]) != 511 && abs(p_pid[p_M1[p_i]]) != 521 && abs(p_pid[p_M1[p_i]]) != 531 && abs(p_pid[p_M1[p_i]]) != 541 && abs(p_pid[p_M1[p_i]]) != 5122) {
                        if (abs(p_pid[p_i]) == 411 || abs(p_pid[p_i]) == 421 || abs(p_pid[p_i]) == 431 || abs(p_pid[p_i]) == 4122) {
                                int dummy_mother = p_i;
                                while (abs(p_pid[dummy_mother]) != 23 && abs(p_pid[dummy_mother]) != 511 && abs(p_pid[dummy_mother]) != 521 && abs(p_pid[dummy_mother]) != 531 && abs(p_pid[dummy_mother]) != 541 && abs(p_pid[dummy_mother]) != 5122) {
                                        dummy_mother = p_M1[dummy_mother];
                                }
                                if (abs(p_pid[dummy_mother]) == 23) {
                                        c_size += 1;
                                }


                                if (c1_idx == 0) { c1_idx = p_i; }
                                if (c1_idx != 0) { c2_idx = p_i; }

                        }
                }

                if (abs(p_pid[p_i]) == 13 && abs(p_pid[p_M1[p_i]]) != 541 && abs(p_pid[p_M1[p_i]]) != 443 && Bc_size == 0) {
                        background_found = 1;

                }
        }

        int two_B = 0, four_B = 0, two_c = 0;
        int B_size = B0_size + Bpm_size + Bs0_size + Lambda_size + Xi_size;
        if (background_found == 0) { continue; }
        //cout << i_en << endl;
        //cout << B_size << "; " << c_size << endl;
        if (background_found == 1 && Bc_size == 0 && B_size == 2 && c_size == 0) { two_B = 1; }
        if (background_found == 1 && Bc_size == 0 && B_size == 2 && c_size == 2) { two_c = 1; }
        if (background_found == 1 && Bc_size == 0 && B_size == 4) { four_B = 1; }

        cout << "=== Event " << i_en << " ===" << endl;
        cout << "Result from previous code: " << endl;
        cout << "2b:" << two_B << "; 2b2c:" << two_c << "; 4b:" << four_B << endl;
        cout << endl;
        //cout << c1_idx << "; " << c2_idx << endl;

        //if (two_c != 1) { continue; }
        */

        cout << "////////////////////////////////////////////////////////////////////" << endl;

        //============= printing true decay ===============
        // to store the indeces of 2b, 4b, 2b2c
        int num_from_Z = 0;
        int levels[100];
        int i_level = 0;
        // find B mesons
        for (int p_i = 0; p_i < p_size; p_i++) {
            int i_mother = p_i;
            // B pid: 5xx, 105xx, 205xx
            int have_B = 0;
            if (abs(int((p_pid[p_i] % 10000) / 100)) == 5) {
                have_B = 1;
            }
            if (have_B == 0) {
                continue;
            }
            int level = 0;
            while (true) {
                i_mother = p_M1[i_mother];
                level += 1;
                if (abs(p_pid[i_mother]) == 23) {
                    break;
                }
            }
            cout << "Found B meson: " << p_pid[p_i] << "\t(product after " << level << " decays from Z)" << endl;
            levels[i_level] = level;
            i_level += 1;
        }

        std::list<double> levels_list(levels, levels + i_level);
        levels_list.sort();
        levels_list.unique();
        int level_arr[levels_list.size()];
        std::copy(levels_list.begin(), levels_list.end(), level_arr);

        // int from_Z = 0;

        int idx_lt[999] = {0};
        int _num_b_ = 0;
        int _num_c_ = 0;
        cout << endl;
        cout << "=== Event " << i_en << " ===" << endl;
        for (int i_size = 0; i_size < levels_list.size(); i_size++) {
            cout << endl;
            for (int p_i = 0; p_i < p_size; p_i++) {
                string chain;
                int i_mother = p_i;
                int from_QCD = 1;
                for (int _ = 0; _ < level_arr[i_size]; _++) {
                    string chain_i = std::to_string(p_pid[i_mother]);
                    i_mother = p_M1[i_mother];
                    chain = chain + chain_i + "\t< ";
                    if (abs(int(p_pid[i_mother] / 100)) != 0) {
                        from_QCD = 0;
                        break;
                    }
                }
                if (from_QCD == 0) {
                    continue;
                }
                chain = chain + "23";

                if (abs(p_pid[i_mother]) != 23) {
                    continue;
                }
                // c pid: 4xx, 104xx, 204xx
                if (abs(int((p_pid[p_i] % 10000) / 100)) != 5 && abs(int((p_pid[p_i] % 10000) / 100)) != 4) {
                    continue;
                }
                if (abs(int((p_pid[p_i] % 10000) / 100)) == 5) {
                    idx_lt[_num_b_ + _num_c_] = p_i;
                    _num_b_ += 1;
                }
                if (abs(int((p_pid[p_i] % 10000) / 100)) == 4) {
                    idx_lt[_num_b_ + _num_c_] = p_i;
                    _num_c_ += 1;
                }

                TLorentzVector meson;
                meson.SetPtEtaPhiE(p_pt[p_i], p_eta[p_i], p_phi[p_i], p_E[p_i]);
                // chain = std::to_string(level_arr[i_size]) + "-th product: (P=" + std::to_string(meson.P()) + ")\t" + chain;
                chain = std::to_string(level_arr[i_size]) + "-th product: (E=" + std::to_string(meson.E()) + ")\t" + chain;

                cout << chain << endl;
            }
        }

        // cout << _num_b_ + _num_c_ << endl;
        Float_t angles[999] = {999};
        Float_t energies[999] = {999};
        int i_angles = 0;
        int i_energies = 0;
        for (int i_qk1 = 0; i_qk1 < _num_b_ + _num_c_; i_qk1++) {
            TLorentzVector meson1;
            meson1.SetPtEtaPhiE(p_pt[idx_lt[i_qk1]], p_eta[idx_lt[i_qk1]], p_phi[idx_lt[i_qk1]], p_E[idx_lt[i_qk1]]);
            for (int i_qk2 = 0; i_qk2 < i_qk1; i_qk2++) {
                TLorentzVector meson2;
                meson2.SetPtEtaPhiE(p_pt[idx_lt[i_qk2]], p_eta[idx_lt[i_qk2]], p_phi[idx_lt[i_qk2]], p_E[idx_lt[i_qk2]]);
                Float_t theta = acos((meson1.Px() * meson2.Px() + meson1.Py() * meson2.Py() + meson1.Pz() * meson2.Pz()) / (meson1.P() * meson2.P()));
                cout << "Angle Bewteen " << p_pid[idx_lt[i_qk1]] << " & " << p_pid[idx_lt[i_qk2]] << ": " << theta * 180.0 / 3.141592654 << " degrees." << endl;
                angles[i_angles] = theta * 180.0 / 3.141592654;
                i_angles += 1;
            }
            energies[i_energies] = meson1.E();
            i_energies += 1;
        }
        sort(angles, angles + i_angles);
        sort(energies, energies + i_energies);

        cout << _num_b_ << "; " << _num_c_ << endl;
        cout << "////////////////////////////////////////////////////////////////////" << endl;

        env_idx_store = i_en;
        num_b_store = _num_b_;
        num_c_store = _num_c_;

        // cout << num_c_store << endl;
        angle1_store = angles[0];
        angle2_store = angles[1];
        E1_store = energies[0];
        E2_store = energies[1];
        // cout << E1_store << "; " << E2_store << "; " << endl;

        int bkg_type;
        if (_num_b_ == 2 && _num_c_ == 2) {
            bkg_type = 1;
        }  // 2b2c
        else if (_num_b_ == 4) {
            bkg_type = 2;
        }  // 4b
        else if (_num_b_ == 2) {
            bkg_type = 2;
        }  // 4b

        // if (false) {
        if (true) {
            // if (bkg_type == 0) {
            // if (bkg_type == 1) {
            // if (bkg_type == 2) {

            Float_t H_EE_0 = 0;
            Float_t H_EE_1 = 0;
            Float_t H_EE_2 = 0;
            Float_t H_EE_3 = 0;
            Float_t H_EE_4 = 0;
            Float_t H_EE_5 = 0;
            Float_t H_EE_6 = 0;
            Float_t H_EE_7 = 0;
            Float_t H_EE_8 = 0;
            Float_t H_EE_9 = 0;
            Float_t H_EE_10 = 0;
            // 1st tracks loop
            for (int i_tr1 = 0; i_tr1 < tr_size; i_tr1++) {
                TLorentzVector cp1;
                cp1.SetPtEtaPhiE(tr_pt[i_tr1], tr_eta[i_tr1], tr_phi[i_tr1], tr_p[i_tr1]);
                Float_t E1 = cp1.E();
                for (int i_tr2 = 0; i_tr2 < tr_size; i_tr2++) {
                    TLorentzVector cp2;
                    cp2.SetPtEtaPhiE(tr_pt[i_tr2], tr_eta[i_tr2], tr_phi[i_tr2], tr_p[i_tr2]);
                    Float_t E2 = cp2.E();
                    Float_t cosine_theta = (cp1.Px() * cp2.Px() + cp1.Py() * cp2.Py() + cp1.Pz() * cp2.Pz()) / (cp1.P() * cp2.P());
                    H_EE_0 += E1 * E2;
                    H_EE_1 += E1 * E2 * cosine_theta;
                    H_EE_2 += E1 * E2 * (0.5 * (3 * pow(cosine_theta, 2) - 1));
                    H_EE_3 += E1 * E2 * (0.5 * (5 * pow(cosine_theta, 3) - 3 * cosine_theta));
                    H_EE_4 += E1 * E2 * ((1.0 / 8.0) * (35 * pow(cosine_theta, 4) - 30 * pow(cosine_theta, 2) + 3));
                    H_EE_5 += E1 * E2 * ((1.0 / 8.0) * (63 * pow(cosine_theta, 5) - 70 * pow(cosine_theta, 3) + 15 * cosine_theta));
                    H_EE_6 += E1 * E2 * ((1.0 / 16.0) * (231 * pow(cosine_theta, 6) - 315 * pow(cosine_theta, 4) + 105 * pow(cosine_theta, 2) - 5));
                    H_EE_7 += E1 * E2 * ((1.0 / 16.0) * (429 * pow(cosine_theta, 7) - 693 * pow(cosine_theta, 5) + 315 * pow(cosine_theta, 3) - 35 * cosine_theta));
                    H_EE_8 += E1 * E2 * ((1.0 / 128.0) * (6435 * pow(cosine_theta, 8) - 12012 * pow(cosine_theta, 6) + 6930 * pow(cosine_theta, 4) - 1260 * pow(cosine_theta, 2) + 35));
                    H_EE_9 += E1 * E2 * ((1.0 / 128.0) * (12155 * pow(cosine_theta, 9) - 25740 * pow(cosine_theta, 7) + 18018 * pow(cosine_theta, 5) - 4620 * pow(cosine_theta, 3) + 315 * cosine_theta));
                    H_EE_10 += E1 * E2 * ((1.0 / 256) * (46189 * pow(cosine_theta, 10) - 109395 * pow(cosine_theta, 8) + 90090 * pow(cosine_theta, 6) - 30030 * pow(cosine_theta, 4) + 3465 * pow(cosine_theta, 2) - 63));
                }
                for (int i_nh2 = 0; i_nh2 < efnh_size; i_nh2++) {
                    Float_t pt2 = efnh_et[i_nh2] * pow(efnh_E[i_nh2] * efnh_E[i_nh2] - m_K0 * m_K0, 0.5) / efnh_E[i_nh2];
                    TLorentzVector nh2;
                    nh2.SetPtEtaPhiE(pt2, efnh_eta[i_nh2], efnh_phi[i_nh2], efnh_E[i_nh2]);
                    Float_t E2 = nh2.E();
                    Float_t cosine_theta = (cp1.Px() * nh2.Px() + cp1.Py() * nh2.Py() + cp1.Pz() * nh2.Pz()) / (cp1.P() * nh2.P());
                    H_EE_0 += E1 * E2;
                    H_EE_1 += E1 * E2 * cosine_theta;
                    H_EE_2 += E1 * E2 * (0.5 * (3 * pow(cosine_theta, 2) - 1));
                    H_EE_3 += E1 * E2 * (0.5 * (5 * pow(cosine_theta, 3) - 3 * cosine_theta));
                    H_EE_4 += E1 * E2 * ((1.0 / 8.0) * (35 * pow(cosine_theta, 4) - 30 * pow(cosine_theta, 2) + 3));
                    H_EE_5 += E1 * E2 * ((1.0 / 8.0) * (63 * pow(cosine_theta, 5) - 70 * pow(cosine_theta, 3) + 15 * cosine_theta));
                    H_EE_6 += E1 * E2 * ((1.0 / 16.0) * (231 * pow(cosine_theta, 6) - 315 * pow(cosine_theta, 4) + 105 * pow(cosine_theta, 2) - 5));
                    H_EE_7 += E1 * E2 * ((1.0 / 16.0) * (429 * pow(cosine_theta, 7) - 693 * pow(cosine_theta, 5) + 315 * pow(cosine_theta, 3) - 35 * cosine_theta));
                    H_EE_8 += E1 * E2 * ((1.0 / 128.0) * (6435 * pow(cosine_theta, 8) - 12012 * pow(cosine_theta, 6) + 6930 * pow(cosine_theta, 4) - 1260 * pow(cosine_theta, 2) + 35));
                    H_EE_9 += E1 * E2 * ((1.0 / 128.0) * (12155 * pow(cosine_theta, 9) - 25740 * pow(cosine_theta, 7) + 18018 * pow(cosine_theta, 5) - 4620 * pow(cosine_theta, 3) + 315 * cosine_theta));
                    H_EE_10 += E1 * E2 * ((1.0 / 256) * (46189 * pow(cosine_theta, 10) - 109395 * pow(cosine_theta, 8) + 90090 * pow(cosine_theta, 6) - 30030 * pow(cosine_theta, 4) + 3465 * pow(cosine_theta, 2) - 63));
                }
                for (int i_ph2 = 0; i_ph2 < efp_size; i_ph2++) {
                    TLorentzVector ph2;
                    ph2.SetPtEtaPhiE(efp_et[i_ph2], efp_eta[i_ph2], efp_phi[i_ph2], efp_E[i_ph2]);
                    Float_t E2 = ph2.E();
                    Float_t cosine_theta = (cp1.Px() * ph2.Px() + cp1.Py() * ph2.Py() + cp1.Pz() * ph2.Pz()) / (cp1.P() * ph2.P());
                    H_EE_0 += E1 * E2;
                    H_EE_1 += E1 * E2 * cosine_theta;
                    H_EE_2 += E1 * E2 * (0.5 * (3 * pow(cosine_theta, 2) - 1));
                    H_EE_3 += E1 * E2 * (0.5 * (5 * pow(cosine_theta, 3) - 3 * cosine_theta));
                    H_EE_4 += E1 * E2 * ((1.0 / 8.0) * (35 * pow(cosine_theta, 4) - 30 * pow(cosine_theta, 2) + 3));
                    H_EE_5 += E1 * E2 * ((1.0 / 8.0) * (63 * pow(cosine_theta, 5) - 70 * pow(cosine_theta, 3) + 15 * cosine_theta));
                    H_EE_6 += E1 * E2 * ((1.0 / 16.0) * (231 * pow(cosine_theta, 6) - 315 * pow(cosine_theta, 4) + 105 * pow(cosine_theta, 2) - 5));
                    H_EE_7 += E1 * E2 * ((1.0 / 16.0) * (429 * pow(cosine_theta, 7) - 693 * pow(cosine_theta, 5) + 315 * pow(cosine_theta, 3) - 35 * cosine_theta));
                    H_EE_8 += E1 * E2 * ((1.0 / 128.0) * (6435 * pow(cosine_theta, 8) - 12012 * pow(cosine_theta, 6) + 6930 * pow(cosine_theta, 4) - 1260 * pow(cosine_theta, 2) + 35));
                    H_EE_9 += E1 * E2 * ((1.0 / 128.0) * (12155 * pow(cosine_theta, 9) - 25740 * pow(cosine_theta, 7) + 18018 * pow(cosine_theta, 5) - 4620 * pow(cosine_theta, 3) + 315 * cosine_theta));
                    H_EE_10 += E1 * E2 * ((1.0 / 256) * (46189 * pow(cosine_theta, 10) - 109395 * pow(cosine_theta, 8) + 90090 * pow(cosine_theta, 6) - 30030 * pow(cosine_theta, 4) + 3465 * pow(cosine_theta, 2) - 63));
                }
            }
            // 1st neutral hadron loop
            for (int i_nh1 = 0; i_nh1 < efnh_size; i_nh1++) {
                Float_t pt1 = efnh_et[i_nh1] * pow(efnh_E[i_nh1] * efnh_E[i_nh1] - m_K0 * m_K0, 0.5) / efnh_E[i_nh1];
                TLorentzVector nh1;
                nh1.SetPtEtaPhiE(pt1, efnh_eta[i_nh1], efnh_phi[i_nh1], efnh_E[i_nh1]);
                Float_t E1 = nh1.E();
                // Float_t E1 = nh1.P();
                for (int i_tr2 = 0; i_tr2 < tr_size; i_tr2++) {
                    TLorentzVector cp2;
                    cp2.SetPtEtaPhiE(tr_pt[i_tr2], tr_eta[i_tr2], tr_phi[i_tr2], tr_p[i_tr2]);
                    Float_t E2 = cp2.E();
                    // Float_t E2 = cp2.P();
                    Float_t cosine_theta = (nh1.Px() * cp2.Px() + nh1.Py() * cp2.Py() + nh1.Pz() * cp2.Pz()) / (nh1.P() * cp2.P());
                    H_EE_0 += E1 * E2;
                    H_EE_1 += E1 * E2 * cosine_theta;
                    H_EE_2 += E1 * E2 * (0.5 * (3 * pow(cosine_theta, 2) - 1));
                    H_EE_3 += E1 * E2 * (0.5 * (5 * pow(cosine_theta, 3) - 3 * cosine_theta));
                    H_EE_4 += E1 * E2 * ((1.0 / 8.0) * (35 * pow(cosine_theta, 4) - 30 * pow(cosine_theta, 2) + 3));
                    H_EE_5 += E1 * E2 * ((1.0 / 8.0) * (63 * pow(cosine_theta, 5) - 70 * pow(cosine_theta, 3) + 15 * cosine_theta));
                    H_EE_6 += E1 * E2 * ((1.0 / 16.0) * (231 * pow(cosine_theta, 6) - 315 * pow(cosine_theta, 4) + 105 * pow(cosine_theta, 2) - 5));
                    H_EE_7 += E1 * E2 * ((1.0 / 16.0) * (429 * pow(cosine_theta, 7) - 693 * pow(cosine_theta, 5) + 315 * pow(cosine_theta, 3) - 35 * cosine_theta));
                    H_EE_8 += E1 * E2 * ((1.0 / 128.0) * (6435 * pow(cosine_theta, 8) - 12012 * pow(cosine_theta, 6) + 6930 * pow(cosine_theta, 4) - 1260 * pow(cosine_theta, 2) + 35));
                    H_EE_9 += E1 * E2 * ((1.0 / 128.0) * (12155 * pow(cosine_theta, 9) - 25740 * pow(cosine_theta, 7) + 18018 * pow(cosine_theta, 5) - 4620 * pow(cosine_theta, 3) + 315 * cosine_theta));
                    H_EE_10 += E1 * E2 * ((1.0 / 256) * (46189 * pow(cosine_theta, 10) - 109395 * pow(cosine_theta, 8) + 90090 * pow(cosine_theta, 6) - 30030 * pow(cosine_theta, 4) + 3465 * pow(cosine_theta, 2) - 63));
                }
                for (int i_nh2 = 0; i_nh2 < efnh_size; i_nh2++) {
                    Float_t pt2 = efnh_et[i_nh2] * pow(efnh_E[i_nh2] * efnh_E[i_nh2] - m_K0 * m_K0, 0.5) / efnh_E[i_nh2];
                    TLorentzVector nh2;
                    nh2.SetPtEtaPhiE(pt2, efnh_eta[i_nh2], efnh_phi[i_nh2], efnh_E[i_nh2]);
                    Float_t E2 = nh2.E();
                    // Float_t E2 = nh2.P();
                    Float_t cosine_theta = (nh1.Px() * nh2.Px() + nh1.Py() * nh2.Py() + nh1.Pz() * nh2.Pz()) / (nh1.P() * nh2.P());
                    H_EE_0 += E1 * E2;
                    H_EE_1 += E1 * E2 * cosine_theta;
                    H_EE_2 += E1 * E2 * (0.5 * (3 * pow(cosine_theta, 2) - 1));
                    H_EE_3 += E1 * E2 * (0.5 * (5 * pow(cosine_theta, 3) - 3 * cosine_theta));
                    H_EE_4 += E1 * E2 * ((1.0 / 8.0) * (35 * pow(cosine_theta, 4) - 30 * pow(cosine_theta, 2) + 3));
                    H_EE_5 += E1 * E2 * ((1.0 / 8.0) * (63 * pow(cosine_theta, 5) - 70 * pow(cosine_theta, 3) + 15 * cosine_theta));
                    H_EE_6 += E1 * E2 * ((1.0 / 16.0) * (231 * pow(cosine_theta, 6) - 315 * pow(cosine_theta, 4) + 105 * pow(cosine_theta, 2) - 5));
                    H_EE_7 += E1 * E2 * ((1.0 / 16.0) * (429 * pow(cosine_theta, 7) - 693 * pow(cosine_theta, 5) + 315 * pow(cosine_theta, 3) - 35 * cosine_theta));
                    H_EE_8 += E1 * E2 * ((1.0 / 128.0) * (6435 * pow(cosine_theta, 8) - 12012 * pow(cosine_theta, 6) + 6930 * pow(cosine_theta, 4) - 1260 * pow(cosine_theta, 2) + 35));
                    H_EE_9 += E1 * E2 * ((1.0 / 128.0) * (12155 * pow(cosine_theta, 9) - 25740 * pow(cosine_theta, 7) + 18018 * pow(cosine_theta, 5) - 4620 * pow(cosine_theta, 3) + 315 * cosine_theta));
                    H_EE_10 += E1 * E2 * ((1.0 / 256) * (46189 * pow(cosine_theta, 10) - 109395 * pow(cosine_theta, 8) + 90090 * pow(cosine_theta, 6) - 30030 * pow(cosine_theta, 4) + 3465 * pow(cosine_theta, 2) - 63));
                }
                for (int i_ph2 = 0; i_ph2 < efp_size; i_ph2++) {
                    TLorentzVector ph2;
                    ph2.SetPtEtaPhiE(efp_et[i_ph2], efp_eta[i_ph2], efp_phi[i_ph2], efp_E[i_ph2]);
                    Float_t E2 = ph2.E();
                    // Float_t E2 = ph2.P();
                    Float_t cosine_theta = (nh1.Px() * ph2.Px() + nh1.Py() * ph2.Py() + nh1.Pz() * ph2.Pz()) / (nh1.P() * ph2.P());
                    H_EE_0 += E1 * E2;
                    H_EE_1 += E1 * E2 * cosine_theta;
                    H_EE_2 += E1 * E2 * (0.5 * (3 * pow(cosine_theta, 2) - 1));
                    H_EE_3 += E1 * E2 * (0.5 * (5 * pow(cosine_theta, 3) - 3 * cosine_theta));
                    H_EE_4 += E1 * E2 * ((1.0 / 8.0) * (35 * pow(cosine_theta, 4) - 30 * pow(cosine_theta, 2) + 3));
                    H_EE_5 += E1 * E2 * ((1.0 / 8.0) * (63 * pow(cosine_theta, 5) - 70 * pow(cosine_theta, 3) + 15 * cosine_theta));
                    H_EE_6 += E1 * E2 * ((1.0 / 16.0) * (231 * pow(cosine_theta, 6) - 315 * pow(cosine_theta, 4) + 105 * pow(cosine_theta, 2) - 5));
                    H_EE_7 += E1 * E2 * ((1.0 / 16.0) * (429 * pow(cosine_theta, 7) - 693 * pow(cosine_theta, 5) + 315 * pow(cosine_theta, 3) - 35 * cosine_theta));
                    H_EE_8 += E1 * E2 * ((1.0 / 128.0) * (6435 * pow(cosine_theta, 8) - 12012 * pow(cosine_theta, 6) + 6930 * pow(cosine_theta, 4) - 1260 * pow(cosine_theta, 2) + 35));
                    H_EE_9 += E1 * E2 * ((1.0 / 128.0) * (12155 * pow(cosine_theta, 9) - 25740 * pow(cosine_theta, 7) + 18018 * pow(cosine_theta, 5) - 4620 * pow(cosine_theta, 3) + 315 * cosine_theta));
                    H_EE_10 += E1 * E2 * ((1.0 / 256) * (46189 * pow(cosine_theta, 10) - 109395 * pow(cosine_theta, 8) + 90090 * pow(cosine_theta, 6) - 30030 * pow(cosine_theta, 4) + 3465 * pow(cosine_theta, 2) - 63));
                }
            }
            // 1st photon loop
            for (int i_ph1 = 0; i_ph1 < efp_size; i_ph1++) {
                TLorentzVector ph1;
                ph1.SetPtEtaPhiE(efp_et[i_ph1], efp_eta[i_ph1], efp_phi[i_ph1], efp_E[i_ph1]);
                Float_t E1 = ph1.E();
                // Float_t E1 = ph1.P();
                for (int i_tr2 = 0; i_tr2 < tr_size; i_tr2++) {
                    TLorentzVector cp2;
                    cp2.SetPtEtaPhiE(tr_pt[i_tr2], tr_eta[i_tr2], tr_phi[i_tr2], tr_p[i_tr2]);
                    Float_t E2 = cp2.E();
                    // Float_t E2 = cp2.P();
                    Float_t cosine_theta = (ph1.Px() * cp2.Px() + ph1.Py() * cp2.Py() + ph1.Pz() * cp2.Pz()) / (ph1.P() * cp2.P());
                    H_EE_0 += E1 * E2;
                    H_EE_1 += E1 * E2 * cosine_theta;
                    H_EE_2 += E1 * E2 * (0.5 * (3 * pow(cosine_theta, 2) - 1));
                    H_EE_3 += E1 * E2 * (0.5 * (5 * pow(cosine_theta, 3) - 3 * cosine_theta));
                    H_EE_4 += E1 * E2 * ((1.0 / 8.0) * (35 * pow(cosine_theta, 4) - 30 * pow(cosine_theta, 2) + 3));
                    H_EE_5 += E1 * E2 * ((1.0 / 8.0) * (63 * pow(cosine_theta, 5) - 70 * pow(cosine_theta, 3) + 15 * cosine_theta));
                    H_EE_6 += E1 * E2 * ((1.0 / 16.0) * (231 * pow(cosine_theta, 6) - 315 * pow(cosine_theta, 4) + 105 * pow(cosine_theta, 2) - 5));
                    H_EE_7 += E1 * E2 * ((1.0 / 16.0) * (429 * pow(cosine_theta, 7) - 693 * pow(cosine_theta, 5) + 315 * pow(cosine_theta, 3) - 35 * cosine_theta));
                    H_EE_8 += E1 * E2 * ((1.0 / 128.0) * (6435 * pow(cosine_theta, 8) - 12012 * pow(cosine_theta, 6) + 6930 * pow(cosine_theta, 4) - 1260 * pow(cosine_theta, 2) + 35));
                    H_EE_9 += E1 * E2 * ((1.0 / 128.0) * (12155 * pow(cosine_theta, 9) - 25740 * pow(cosine_theta, 7) + 18018 * pow(cosine_theta, 5) - 4620 * pow(cosine_theta, 3) + 315 * cosine_theta));
                    H_EE_10 += E1 * E2 * ((1.0 / 256) * (46189 * pow(cosine_theta, 10) - 109395 * pow(cosine_theta, 8) + 90090 * pow(cosine_theta, 6) - 30030 * pow(cosine_theta, 4) + 3465 * pow(cosine_theta, 2) - 63));
                }
                for (int i_nh2 = 0; i_nh2 < efnh_size; i_nh2++) {
                    Float_t pt2 = efnh_et[i_nh2] * pow(efnh_E[i_nh2] * efnh_E[i_nh2] - m_K0 * m_K0, 0.5) / efnh_E[i_nh2];
                    TLorentzVector nh2;
                    nh2.SetPtEtaPhiE(pt2, efnh_eta[i_nh2], efnh_phi[i_nh2], efnh_E[i_nh2]);
                    Float_t E2 = nh2.E();
                    // Float_t E2 = nh2.P();
                    Float_t cosine_theta = (ph1.Px() * nh2.Px() + ph1.Py() * nh2.Py() + ph1.Pz() * nh2.Pz()) / (ph1.P() * nh2.P());
                    H_EE_0 += E1 * E2;
                    H_EE_1 += E1 * E2 * cosine_theta;
                    H_EE_2 += E1 * E2 * (0.5 * (3 * pow(cosine_theta, 2) - 1));
                    H_EE_3 += E1 * E2 * (0.5 * (5 * pow(cosine_theta, 3) - 3 * cosine_theta));
                    H_EE_4 += E1 * E2 * ((1.0 / 8.0) * (35 * pow(cosine_theta, 4) - 30 * pow(cosine_theta, 2) + 3));
                    H_EE_5 += E1 * E2 * ((1.0 / 8.0) * (63 * pow(cosine_theta, 5) - 70 * pow(cosine_theta, 3) + 15 * cosine_theta));
                    H_EE_6 += E1 * E2 * ((1.0 / 16.0) * (231 * pow(cosine_theta, 6) - 315 * pow(cosine_theta, 4) + 105 * pow(cosine_theta, 2) - 5));
                    H_EE_7 += E1 * E2 * ((1.0 / 16.0) * (429 * pow(cosine_theta, 7) - 693 * pow(cosine_theta, 5) + 315 * pow(cosine_theta, 3) - 35 * cosine_theta));
                    H_EE_8 += E1 * E2 * ((1.0 / 128.0) * (6435 * pow(cosine_theta, 8) - 12012 * pow(cosine_theta, 6) + 6930 * pow(cosine_theta, 4) - 1260 * pow(cosine_theta, 2) + 35));
                    H_EE_9 += E1 * E2 * ((1.0 / 128.0) * (12155 * pow(cosine_theta, 9) - 25740 * pow(cosine_theta, 7) + 18018 * pow(cosine_theta, 5) - 4620 * pow(cosine_theta, 3) + 315 * cosine_theta));
                    H_EE_10 += E1 * E2 * ((1.0 / 256) * (46189 * pow(cosine_theta, 10) - 109395 * pow(cosine_theta, 8) + 90090 * pow(cosine_theta, 6) - 30030 * pow(cosine_theta, 4) + 3465 * pow(cosine_theta, 2) - 63));
                }
                for (int i_ph2 = 0; i_ph2 < efp_size; i_ph2++) {
                    TLorentzVector ph2;
                    ph2.SetPtEtaPhiE(efp_et[i_ph2], efp_eta[i_ph2], efp_phi[i_ph2], efp_E[i_ph2]);
                    Float_t E2 = ph2.E();
                    // Float_t E2 = ph2.P();
                    Float_t cosine_theta = (ph1.Px() * ph2.Px() + ph1.Py() * ph2.Py() + ph1.Pz() * ph2.Pz()) / (ph1.P() * ph2.P());
                    H_EE_0 += E1 * E2;
                    H_EE_1 += E1 * E2 * cosine_theta;
                    H_EE_2 += E1 * E2 * (0.5 * (3 * pow(cosine_theta, 2) - 1));
                    H_EE_3 += E1 * E2 * (0.5 * (5 * pow(cosine_theta, 3) - 3 * cosine_theta));
                    H_EE_4 += E1 * E2 * ((1.0 / 8.0) * (35 * pow(cosine_theta, 4) - 30 * pow(cosine_theta, 2) + 3));
                    H_EE_5 += E1 * E2 * ((1.0 / 8.0) * (63 * pow(cosine_theta, 5) - 70 * pow(cosine_theta, 3) + 15 * cosine_theta));
                    H_EE_6 += E1 * E2 * ((1.0 / 16.0) * (231 * pow(cosine_theta, 6) - 315 * pow(cosine_theta, 4) + 105 * pow(cosine_theta, 2) - 5));
                    H_EE_7 += E1 * E2 * ((1.0 / 16.0) * (429 * pow(cosine_theta, 7) - 693 * pow(cosine_theta, 5) + 315 * pow(cosine_theta, 3) - 35 * cosine_theta));
                    H_EE_8 += E1 * E2 * ((1.0 / 128.0) * (6435 * pow(cosine_theta, 8) - 12012 * pow(cosine_theta, 6) + 6930 * pow(cosine_theta, 4) - 1260 * pow(cosine_theta, 2) + 35));
                    H_EE_9 += E1 * E2 * ((1.0 / 128.0) * (12155 * pow(cosine_theta, 9) - 25740 * pow(cosine_theta, 7) + 18018 * pow(cosine_theta, 5) - 4620 * pow(cosine_theta, 3) + 315 * cosine_theta));
                    H_EE_10 += E1 * E2 * ((1.0 / 256) * (46189 * pow(cosine_theta, 10) - 109395 * pow(cosine_theta, 8) + 90090 * pow(cosine_theta, 6) - 30030 * pow(cosine_theta, 4) + 3465 * pow(cosine_theta, 2) - 63));
                }
            }

            // cout << H_EE_1 / pow(91.0, 2) << ", " << H_EE_2 / pow(91.0, 2) << endl;
            // cout << endl;
            // hist_H_EE_1->Fill(H_EE_1 / pow(91.2, 2));
            // hist_H_EE_2->Fill(H_EE_2 / pow(91.2, 2));
            // hist_H_EE_3->Fill(H_EE_3 / pow(91.2, 2));
            // hist_H_EE_4->Fill(H_EE_4 / pow(91.2, 2));

            H_EE0_store = H_EE_0 / pow(91.2, 2);
            H_EE1_store = H_EE_1 / pow(91.2, 2);
            H_EE2_store = H_EE_2 / pow(91.2, 2);
            H_EE3_store = H_EE_3 / pow(91.2, 2);
            H_EE4_store = H_EE_4 / pow(91.2, 2);
            H_EE5_store = H_EE_5 / pow(91.2, 2);
            H_EE6_store = H_EE_6 / pow(91.2, 2);
            H_EE7_store = H_EE_7 / pow(91.2, 2);
            H_EE8_store = H_EE_8 / pow(91.2, 2);
            H_EE9_store = H_EE_9 / pow(91.2, 2);
            H_EE10_store = H_EE_10 / pow(91.2, 2);

            // tr.Fill();
            tr2.Fill();
        }

        // cout << endl;
        // cout << "-----------------------------------------------------------" << endl;
    }

    // c->cd(1);
    // hist_H_EE_1->Draw();
    // c->cd(2);
    // hist_H_EE_2->Draw();
    // c->cd(3);
    // hist_H_EE_3->Draw();
    // c->cd(4);
    // hist_H_EE_4->Draw();

    // tr.Write();
    // fea.Close();

    // tr2.Write();
    // fea2.Close();
}
