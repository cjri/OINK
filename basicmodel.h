
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>

struct run_params {
    // Parameters for this simulation run
    double incubation_a; // Parameter for Weibull distribution modelling 
    double incubation_b;
    double infection_a;
    double infection_b;
    int max_infections; // Limit on how many infections to generate
    int infection_length;  
    int add_limit; //Maximum number of infections per generation.
    double probability_detect;   // P(see doctor and get tested if sick)
    double probability_first_detect; // Probability for detecting first case may be different e.g. if known
    int time_symptom_onset_to_detect;
    int replicas; // Number of simulations to perform for each value of R0
    double max_R0; //
    std::vector<double> R0_vals; 
    int max_simulation_time; // Maximum length of simulation 
    int more_stats; //Additional files outputted
    int verb; // Whether to be verbose
    int seed; // Fixed RNG seed
    std::string output_prefix;
    std::string input_prefix;

};

struct patient {
    int time_infected; // Time of infection (relative to index case)
    int time_symptom_onset; //Symptom onset (i.e. non-latent)
    unsigned int detected; // Whether case was detected
    int time_reported; // Time esult reported
};

struct outbreak {
    std::vector<patient> individuals;
    int time_first_detect; // Time of first detection, relative to start of simulation
    int last_time_simulated; // Time at which simulation stopped.
};

struct detect {
    int day;
    int cases;
};

struct output { // Structure to combine output over simulations for each R0 and timepoint
    int tested;
    int accepted;
    int dead;
    std::vector<int> current_size;
    std::vector<int> origin_time;
};
