#include "simulation.hpp"

#include "TAxis.h"
#include "TMath.h"

using namespace std;

Simulation::Simulation(string name, Double_t freq, Double_t pn, Double_t pe,
                       Double_t c, Double_t temperature,
                       Double_t t1n, Double_t t1e)
    : name(name), t1n(t1n), t1e(t1e), c(c), phi(0.0), pn_raw(pn), pn(pn), pe(pe), dose(0.0) {
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

void Simulation::beam_on(Double_t current) {
    beam_current = current;
}

void Simulation::beam_off() {
    beam_current = 0.0;
}

void Simulation::run_until(Double_t t_final) {
    while (t < t_final) {
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
    // Parameters for temperature change (exponential growth/decay)
    // TEMP_SS = steady-state temperature
    // K_TEMP = rate of exponential increase
    const auto K_TEMP = 1.0;
    const auto TEMP_SS = 1.0 + beam_current / 100.0;

    // Increase phi according to some exponential growth when the beam is on
    const auto K_PHI = beam_current / 1e7;
    const auto PHI_SS = 0.001;

    for (int i = 0; i < N_ITER; i++) {
        // Calculate constants (for convenience)
        const auto a_const = -t1e / t1n - (c / 2) * (alpha + beta) - phi;
        const auto b_const = (c / 2) * (alpha - beta);
        const auto c_const = (alpha - beta) / 2;
        const auto d_const = -1 - (alpha + beta) / 2;

        // Calculate rates
        const auto pn_prime = (a_const * pn_raw + b_const * pe) / t1e;
        const auto pe_prime = (c_const * pn_raw + d_const * pe + pe0) / t1e;

        // Update pn and pe (Euler's method)
        pn_raw += pn_prime * TIME_STEP / N_ITER;
        pe += pe_prime * TIME_STEP / N_ITER;
        // Update temperature and phi
        set_temperature(temperature + (TIME_STEP / N_ITER) * K_TEMP * (TEMP_SS - temperature));
        phi += (TIME_STEP / N_ITER) * K_PHI * (PHI_SS - phi);
    }

    // Update "noisy pn"
    pn = pn_noisy();

    // Add dose
    // Must convert beam_current (in nA) to electrons/second
    dose += (beam_current * 1e-9 / ELEM_CHARGE) * TIME_STEP;

    t += TIME_STEP;
}

Double_t Simulation::pn_noisy() {
    const Double_t noise = RANDOMNESS * (0.5 - rng->Rndm());
    return pn_raw * (1 + noise);
}
