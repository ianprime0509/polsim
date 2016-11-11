#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP

#include "TROOT.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TTree.h"

#include <memory>
#include <string>

/* The following class represents a single solid polarized target experiment
 * in its entirety (i.e. it takes into account all relevant parameters).
 * The simulation stores data internally, which can be written to a ROOT
 * file or graphed to a TCanvas using the draw method.
 *
 * Some notes on the physics involved:
 * BEAM: the beam is activated by using the beam_on method, which takes
 * as its argument a beam current in nA. The beam will both heat
 * the material and change its paramagnetic structure; currently,
 * the latter effect is modeled by increasing the "extraneous relaxation rate"
 * phi according to a certain exponential growth (in reality, it should change only
 * the parameter C, but the underlying physical model does not produce the correct
 * results under this change).
 *
 * Units used in the below quantities:
 * Polarization: as a number in the range [-1, 1]
 * Time: in seconds
 * Temperature: in K
 * Field: in T
 * Dose: in electrons / cm^3
 * T1N, T1E: in seconds
 * Beam current: in nA
 */
class Simulation {
    // Physical parameters
    static constexpr Double_t ELEM_CHARGE = 1.602176662e-19; // in C
    // Parameters for alpha/beta vs frequency fit
    static constexpr Double_t FIT_A = 0.545266;
    static constexpr Double_t FIT_S = 0.088415;
    // The m1 and m2 parameters are the centers of the beta/alpha
    // distributions, respectively. The actual center changes
    // as dose is added, according to:
    // m1 = (FIT_M1_BASE - FIT_M1_COEFF) + FIT_M1_COEFF * exp(FIT_M1_RATE * dose)
    // and similarly for m2.
    // The justification for this model can be found in the 2016 presentation,
    // which itself got these parameters from the SANE frequency vs dose data.
    // The FIT_M1_BASE and FIT_M2_BASE parameters are the centers of the distributions
    // when polarizing with the beam off (in GHz).
    static constexpr Double_t FIT_M1_BASE = 140.286;
    static constexpr Double_t FIT_M1_COEFF = 0.045;
    static constexpr Double_t FIT_M1_RATE = -0.38e-15;
    static constexpr Double_t FIT_M2_BASE = 140.468;
    static constexpr Double_t FIT_M2_COEFF = -0.065;
    static constexpr Double_t FIT_M2_RATE = -3.8e-15;
    // Simulation parameters
    static constexpr Double_t TIME_STEP = 1;
    static constexpr int N_ITER = 1000;
    static constexpr Double_t RANDOMNESS = 0.02;
    // Not sure what kind of a parameter this is
    // Whenever the beam is on, C is increased according to
    // (delta)C = IRRADIATION_FACTOR * beam_current * (delta)t
    static constexpr Double_t IRRADIATION_FACTOR = 1e-10;

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
    // Internal state variables
    bool in_anneal;

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
    // Anneals the material for the given time at the given temperature
    void anneal(Double_t time, Double_t temp=70.0);
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
    // Calculates the alpha and beta parameters
    void calc_transition_rates();
    // Returns the value of pn with some noise applied (thermal fluctuations)
    Double_t pn_noisy();
};

#endif
