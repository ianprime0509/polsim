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
    system_temperature = temperature;
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
    freq = f;
    // Make sure to calculate alpha and beta parameters
    calc_transition_rates();
}

void Simulation::set_system_temperature(Double_t temp) {
    system_temperature = temp;
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

void Simulation::anneal(Double_t time, Double_t temp) {
    // Reset phi (i.e. remove negative effects of irradiation)
    phi = 0.0;
    // Maybe change t1n?
    t1n *= 0.8;

    auto temp_tmp = system_temperature;
    set_system_temperature(temp);
    run_until(t + time);
    set_system_temperature(temp_tmp);
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

void Simulation::set_temperature(Double_t temp) {
    pe0 = -TMath::TanH(2 / temp);
    temperature = temp;
}

void Simulation::time_step() {
    // Parameters for temperature change (exponential growth/decay)
    // TEMP_SS = steady-state temperature
    // K_TEMP = rate of exponential increase
    // If we're annealing, we shouldn't allow the temperature to change
    // (assume anneals occur at constant temperature)
    const auto K_TEMP = 0.01;
    const auto TEMP_SS = system_temperature + beam_current / 100.0;

    // Increase phi according to some exponential growth when the beam is on
    // Parameters are similar to those for temperature change
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

        // Update C and dose
        c += IRRADIATION_FACTOR * beam_current * (TIME_STEP / N_ITER);
        // Must convert beam_current (in nA) to electrons/second
        dose += (beam_current * 1e-9 / ELEM_CHARGE) * TIME_STEP / N_ITER;
        // Calculate new transition rates (alpha and beta)
        calc_transition_rates();

        // Update time
        t += TIME_STEP / N_ITER;
    }

    // Update "noisy pn"
    pn = pn_noisy();
}

void Simulation::calc_transition_rates() {
    // Calculate distribution parameters (the means m1 and m2 are particularly important)
    const auto FIT_M1 = (FIT_M1_BASE - FIT_M1_COEFF) + FIT_M1_COEFF * TMath::Exp(FIT_M1_RATE * dose);
    const auto FIT_M2 = (FIT_M2_BASE - FIT_M2_COEFF) + FIT_M2_COEFF * TMath::Exp(FIT_M2_RATE * dose);
    const auto scale = FIT_A / (TMath::Sqrt(2 * TMath::Pi()) * FIT_S);
    const auto diff1 = freq - FIT_M1;
    const auto diff2 = freq - FIT_M2;
    const auto exp1 = TMath::Exp(-diff1 * diff1 / (2 * FIT_S * FIT_S));
    const auto exp2 = TMath::Exp(-diff2 * diff2 / (2 * FIT_S * FIT_S));

    alpha = scale * exp2;
    beta = scale * exp1;
}

Double_t Simulation::pn_noisy() {
    // Account for both types of noise
    const Double_t thermal_noise = THERMAL_RANDOMNESS * (0.5 - rng->Rndm());
    const Double_t uniform_noise = BASE_RANDOMNESS * (0.5 - rng->Rndm());
    return pn_raw * (1 + thermal_noise) + uniform_noise;
}
