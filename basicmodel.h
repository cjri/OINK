using namespace std;
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
    double r0;
    double incubation_a;
    double incubation_b;
    double infection_a;
    double infection_b;
    int max_infections; //Limit on how many infections to generate
    int infection_length;
    int add_limit; //Maximum number of infections per generation.
    double detect;   //P(see doctor and get tested if sick)
    double first_detect; //Probability for detecting first case may be different e.g. if known
    double symptom_to_detect;
    int replicas;
    int run_fast;
    int resolution; //Measured in days - how accurately do we know when things happened?
    string species; //Specify virus
    int test; //Testing the code
    int more_stats; //Additional files outputted
    int verb;
};

//Time relative to first detection?  Could do this retrospectively i.e. once the outbreak is fully simulated?
//What is an efficient way to filter out simulations with too many cases detected?
struct pat {
    int time_i; //Infection
    int time_s; //Symptom onset (i.e. non-latent)
    int infected_by;
    vector<int> infects;
    int detected;
    int see_gp;
    int time_r; //Result reported
};

struct outbreak {
    vector<pat> indiv;
    int first_detect; //Time of first detection.  Update after each generation if not yet found.  Once found shift all times.
    int origin_time;
    int total_detections;
    int last_time_completed; //Starts at zero.  After each generation becomes the time before which nothing else can happen.
    //Idea here: If there are too many detections by last_time_completed, can stop the simulation
    //Implication: Need to have data describing number of detections by what time.
};

struct detect {
    int day;
    int cases;
};

struct output { //Specific to R0
    int tested;
    int accepted;
    int dead;
    vector<int> current_size;
    vector<int> origin_time;
};
