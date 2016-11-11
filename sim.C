// TODO:
// Temperature change is a little arbitrary; what is the actual thermal behavior?
// (e.g. find a good value for K_TEMP)

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

    auto sim1 = new Simulation("sim1", 140.2);
    sim1->run_until(10000);
    sim1->beam_on();
    sim1->run_until(12000);
    sim1->beam_off();
    sim1->run_until(12100);
    sim1->beam_on();
    sim1->run_until(15000);
    sim1->beam_off();
    sim1->run_until(17000);
    sim1->anneal(500);
    sim1->run_until(25000);
    sim1->beam_on();
    sim1->run_until(30000);

    auto sim2 = new Simulation("sim2", 140.2);
    sim2->set_temperature(70.0);
    sim2->run_until(30000);
    sim2->write_data();

    auto canvas = new TCanvas("c1", "Polarization vs time");

    sim1->draw();
    sim2->draw("C same", kRed);

    return 0;
}
