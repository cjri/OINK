
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
    double incubation_a; // Parameter for Weibull distribution modelling incubation length
    double incubation_b; // Parameter for Weibull distribution modelling incubation length
    double infection_a; // Parameter for Weibull distribution modelling time to infection
    double infection_b; // Parameter for Weibull distribution modelling time to infection
    int infection_length_min;  // Number of days a case is considered to be infected. Modelled as uniform on [min, max]
    int infection_length_max;  // Number of days a case is considered to be infected. Modelled as uniform on [min, max]
  
    int add_limit; //Maximum number of infections per generation.
    double probability_detect;   // P(see doctor and get tested if sick)
    double probability_detect_after_first; // P(see doctor and get tested if sick) after first case detected (must be bigger than probability_detect)
    double probability_first_detect; // Probability for detecting first case. By default zero
  
    int time_symptom_onset_to_detect_min; // Number of days from symptom onset to detection. Modelled as uniform on [min, max]
    int time_symptom_onset_to_detect_max; // Number of days from symptom onset to detection. Modelled as uniform on [min, max]
    int replicas; // Number of simulations to perform for each value of R0
    double max_R0; // Simulate R0 values from 0.1 up to max_R0
    std::vector<double> R0_vals; // Vector of R0 valies
    int max_simulation_time; // Maximum length of simulation 
    int more_stats; //Additional files outputted
    int verb; // Whether to be verbose
    int seed; // Fixed RNG seed
    std::string output_prefix; // Prefix for output files
    std::string input_prefix; // Prefix for input data files

};

struct patient {
    int time_infected; // Time of infection (relative to index case)
    int time_symptom_onset; // Time of symptom onset 
    unsigned int detected; // Whether case was detected (as first detected case)
    unsigned int detected_after; // Whether case would be detected if higher surveillance was present
    int infection_length; // Length of infection for that patient.
                          // Patient is infected over [time_symptom_onset, time_symptom_onset+infection_length-1]
    int time_reported; // Time result reported
};

struct outbreak {
    std::vector<patient> individuals; // List storing cases
    int time_first_detect; // Time of first detection, relative to start of simulation
    int time_first_detect_symptom_onset; // Time of symptom onset for the first detected case
    int last_time_simulated; // Time at which simulation stopped.
};

struct detect {
    int day; // Date of detection (relative to first detected case)
    int cases; // Number of cases detected on that day
};

struct output { // Structure to combine output over simulations for each R0 and observation timepoint
    int tested; // Number of tested simulations
    int accepted; // Number of accepted simulations
    int dead; // Number of accepted simulations where there are no actively infected individuals at the observation timepoint
    std::vector<int> current_size; // Total number of currently infected individuals at the observation timepoint
    std::vector<int> origin_time; // Times at which the index case became infected (days prior to the first detection)
};
