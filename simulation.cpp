#include "simulation.hpp"

#include "TAxis.h"
#include "TMath.h"

using namespace std;

Simulation::Simulation(string name, Double_t freq, Double_t pn, Double_t pe,
                       Double_t c, Double_t temperature,
                       Double_t t1n, Double_t t1e)
    : name(name), t1n(t1n), t1e(t1e), c(c), pn_raw(pn), pn(pn), pe(pe), dose(0.0) {
    set_freq(freq);
    set_temperature(temperature);
    beam_off();
    // Set up internal TTree
    tree = unique_ptr<TTree>(new TTree(name.c_str(), "Polarization Data"));
    tree->Branch("time", &this->t, "time/D");
    tree->Branch("pn", &this->pn, "pn/D");
    tree->Branch("pe", &this->pe, "pe/D");
    tree->Branch("frequency", &this->freq, "freq/D");
    tree->Branch("c", &this->c, "c/D");
    tree->Branch("temperature", &this->temperature, "temperature/D");
    tree->Branch("dose", &this->dose, "dose/D");
    // Initialize random number generator
    rng = unique_ptr<TRandom3>(new TRandom3());
}

void Simulation::set_freq(Double_t f) {
    // Make sure to calculate alpha and beta parameters
    const Double_t scale = FIT_A / (TMath::Sqrt(2 * TMath::Pi()) * FIT_S);
    const Double_t diff1 = f - FIT_M1;
    const Double_t diff2 = f - FIT_M2;
    const Double_t exp1 = TMath::Exp(-diff1 * diff1 / (2 * FIT_S * FIT_S));
    const Double_t exp2 = TMath::Exp(-diff2 * diff2 / (2 * FIT_S * FIT_S));

    freq = f;
    alpha = scale * exp2;
    beta = scale * exp1;
}

void Simulation::set_temperature(Double_t temp) {
    pe0 = -TMath::TanH(2 / temp);
    temperature = temp;
}

void Simulation::beam_on() {
    is_beam = true;
}

void Simulation::beam_off() {
    is_beam = false;
}

void Simulation::run_until(Double_t t_final) {
    while (t < t_final) {
        // Update temperature
        // These parameters are different when the beam is on
        // versus when it is off
        // T_SS = steady-state temperature
        // K = rate of exponential increase
        Double_t K, T_SS;
        if (is_beam) {
            K = 100.0;
            T_SS = 5.0;
        } else {
            K = 10.0;
            T_SS = 1;
        }

        auto temp_tmp = temperature;
        // Iterate Euler's method on this temporary temperature
        // This gives better accuracy
        for (auto i = 0; i < N_ITER; i++) {
            temp_tmp += (TIME_STEP / N_ITER) * K * (T_SS - temp_tmp);
        }
        set_temperature(temp_tmp);

        time_step();
        tree->Fill();
    }
}

void Simulation::write_data() {
    tree->Write();
}

void Simulation::draw(const char *options, const Color_t color) {
    Int_t n = tree->GetEntries();
    graph = unique_ptr<TGraph>(new TGraph(n));
    graph->SetTitle((name + " (pn vs t)").c_str());
    graph->GetXaxis()->SetTitle("Time");
    graph->GetYaxis()->SetTitle("Polarization");
    graph->GetXaxis()->Print();
    graph->SetLineColor(color);
    // Save backups since we'll be reading from tree
    auto t_tmp = t;
    auto pn_tmp = pn;
    for (auto i = 0; i < n; i++) {
        tree->GetEntry(i);
        graph->SetPoint(i, t, pn);
    }
    t = t_tmp;
    pn = pn_tmp;

    // Draw graph
    graph->Draw(options);
}

void Simulation::time_step() {
    // Calculate constants
    const auto a_const = -t1e / t1n - (c / 2) * (alpha + beta);
    const auto b_const = (c / 2) * (alpha - beta);
    const auto c_const = (alpha - beta) / 2;
    const auto d_const = -1 - (alpha + beta) / 2;
    for (int i = 0; i < N_ITER; i++) {
        const auto pn_prime = (a_const * pn_raw + b_const * pe) / t1e;
        const auto pe_prime = (c_const * pn_raw + d_const * pe + pe0) / t1e;

        // Update pn and pe (Euler's method)
        pn_raw += pn_prime * TIME_STEP / N_ITER;
        pe += pe_prime * TIME_STEP / N_ITER;
    }

    // Update "noisy pn"
    pn = pn_noisy();

    // Add dose
    if (is_beam) {
        dose += 0.01 * TIME_STEP;
    }

    t += TIME_STEP;
}

Double_t Simulation::pn_noisy() {
    const Double_t noise = RANDOMNESS * (0.5 - rng->Rndm());
    return pn_raw * (1 + noise);
}
