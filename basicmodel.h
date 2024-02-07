
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
    double r0;
    double incubation_a; // Parameter for Weibull distribution modelling 
    double incubation_b;
    double infection_a;
    double infection_b;
    int max_infections; //Limit on how many infections to generate
    int infection_length;
    int add_limit; //Maximum number of infections per generation.
    double probability_detect;   // P(see doctor and get tested if sick)
    double probability_first_detect; // Probability for detecting first case may be different e.g. if known
    double time_symptom_onset_to_detect;
    int replicas;
    int more_stats; //Additional files outputted
    int verb;
    int seed; // Fixed RNG seed
};

struct patient {
    double time_infected; //Infection
    double time_symptom_onset; //Symptom onset (i.e. non-latent)
    int infected_by;
    unsigned int detected;
    int time_reported; //Result reported
};

struct outbreak {
    std::vector<patient> individuals;
    int time_first_detect; //Time of first detection.  Update after each generation if not yet found.  Once found shift all times.
    int origin_time;
    int total_detections;
    int last_time_completed; //Starts at zero.  After each generation becomes the time before which nothing else can happen. (Unclear?)
};

struct detect {
    int day;
    int cases;
};

struct output { //Specific to R0
    run_params params;
    int tested;
    int accepted;
    int dead;
    std::vector<int> current_size;
    std::vector<int> origin_time;
};
