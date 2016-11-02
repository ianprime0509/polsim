// TODO:
// Implement annealing

#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TTree.h"

#include <iostream>
#include <memory>

#include "simulation.hpp"

// TODO: figure out how to actually use ROOT in a normal
// C++ way (linking, anyone?)
#include "simulation.cpp"

using namespace std;

// ROOT entry point
int sim() {
    string filename = "data.root";
    // cout << "Enter filename: " << flush;
    // getline(cin, filename);
    // filename += ".root";
    
    auto output_file = unique_ptr<TFile>(new TFile(filename.c_str(), "recreate"));
    if (output_file->IsOpen()) {
        cout << "Successfully opened file: " << filename << endl;
    } else {
        cerr << "Could not open file: " << filename << endl;
        return 1;
    }

    auto sim0 = new Simulation("sim0", 140.15);
    sim0->run_until(5000);
    sim0->write_data();
    
    auto sim = new Simulation("sim1", 140.15);
    sim->run_until(1000);
    sim->beam_on();
    sim->run_until(1200);
    sim->beam_off();
    sim->run_until(1250);
    sim->beam_on();
    sim->run_until(5000);
    sim->write_data();

    auto sim2 = new Simulation("sim2", 140.15);
    sim2->run_until(1000);
    sim2->beam_on();
    sim2->run_until(5000);
    sim2->write_data();

    auto sim3 = new Simulation("sim3", 140.15);
    sim3->beam_on();
    sim3->run_until(5000);
    sim3->beam_off();
    sim3->write_data();

    auto sim4 = new Simulation("sim4", 140.15);
    sim4->run_until(1000);
    sim4->beam_on();
    sim4->run_until(1200);
    sim4->beam_off();
    sim4->run_until(5000);
    sim4->write_data();

    auto sim5 = new Simulation("sim5", 140.15);
    sim5->beam_on();
    sim5->run_until(3000);
    sim5->beam_off();
    sim5->run_until(5000);

    auto canvas = new TCanvas("c1", "Polarization vs time");
    sim0->draw("AC", kOrange);
    sim->draw("C same");
    sim2->draw("C same", kRed);
    sim3->draw("C same", kBlue);
    sim4->draw("C same", kGreen);
    sim5->draw("C same", kViolet);

    return 0;
}
