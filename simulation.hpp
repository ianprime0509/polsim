#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP

#include "TROOT.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TTree.h"

#include <memory>
#include <string>

class Simulation {
    // Physical parameters
    static constexpr Double_t ELEM_CHARGE = 1.602176662e-19; // in C
    // Parameters for alpha/beta vs frequency fit
    static constexpr Double_t FIT_A = 0.545266;
    static constexpr Double_t FIT_S = 0.088415;
    static constexpr Double_t FIT_M1 = 140.286;
    static constexpr Double_t FIT_M2 = 140.468;
    // Simulation parameters
    static constexpr Double_t TIME_STEP = 1;
    static constexpr int N_ITER = 1000;
    static constexpr Double_t RANDOMNESS = 0.02;

    // Name of the current simulation (for data storage and graphing)
    std::string name;

    // Physical constants
    Double_t t1n, t1e;
    // "External" physical parameters
    Double_t t, freq, temperature;
    // Internal physical parameters
    Double_t alpha, beta, c, pe0, phi;
    // Polarization values
    // pn_raw is the "raw polarization" (without noise)
    Double_t pn_raw, pn, pe;
    // Dose
    Double_t dose, beam_current;

    // TTree (for storing data)
    std::unique_ptr<TTree> tree;
    // TGraph (for plotting the data)
    std::unique_ptr<TGraph> graph;

    // Random number generator
    std::unique_ptr<TRandom3> rng;

public:
    Simulation(std::string name, Double_t freq, Double_t pn = 0.0, Double_t pe = -1.0,
               Double_t c = 0.000136073, Double_t temperature = 1.0,
               Double_t t1n = 25*60.0, Double_t t1e = 0.03);

    void set_freq(Double_t freq);
    void set_temperature(Double_t temperature);
    // Turns the beam on (current measured in nA)
    void beam_on(Double_t current=100.0);
    // Turns the beam off
    void beam_off();
    // Runs the simulation and stops when t >= t_final
    void run_until(Double_t t_final);
    // Writes data to file
    // Why is the file not specified? Because ROOT has a single
    // global file which all writes are directed to!
    // WOW! (who came up with this abomination)
    void write_data();
    // Draws a graph using the given options (as with everything in ROOT,
    // the canvas that it draws to is the global one...
    void draw(const char *options="AC", const Color_t color=kBlack);

private:
    void time_step();
    // Returns the value of pn with some noise applied (thermal fluctuations)
    Double_t pn_noisy();
};

#endif
