#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "basicmodel.h"
#include "io.h"
#include "utilities.h"

int main(int argc, const char **argv)
{

    run_params p;
    GetOptions(p, argc, argv); // Fills run_params structure p with the options for this run

    // Initialize RNG to a fixed seed
    int seed = p.seed;
    gsl_rng_env_setup();
    gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rgen, seed);

    // Import detected infection file
    int n_detections = 0;
    std::vector<detect> detections;
    // Read Data/detections.dat
    ImportDetections(p, n_detections, detections); 
    // detections is a list of detections, with day =integer day index, cases=number of cases on that day.
    // n_detections is the sum of the number of cases.
    
    // We also want to know what time it is, or at least a std::vector of times at which to model what we know.
    int min_time;
    int max_time;
    std::vector<int> timepoints;
    // Read Data/Time_points.dat
    ImportTimePoints(p, min_time, max_time, timepoints);
    // min_time = minimum of timepoints
    // max_time = maximum of timepoints 
    // timepoints = list of integer timepoints (ordered as in file)

    // Generate array of output values by R0 from 0.1 to 4.0
    // For each of these consider each time point.
    std::vector<std::vector<output>> results;
    ConstructResults(p, timepoints, results);
    // For each timepoint, make a list of p.R0_vals.size() zero-initialized output objects
    for (unsigned long r0val=0; r0val<p.R0_vals.size(); r0val++) {
        std::cout << "Generating simulations R0= " << p.R0_vals[r0val] << " "
                  << "\n";
        double r0 = p.R0_vals[r0val];
        for (int calc = 0; calc < p.replicas; calc++)
        { // Generate a number of simulated outbreaks at given R0
            // Construct outbreak shape
            outbreak o;
            SetupOutbreak(o); // Initialize outbreak

            // Set up permutation parameters
            std::vector<int> t_detects_relative; // Detection times ?
            std::vector<int> number_new_symptomatic; // Population sizes ?

            RunSimulationTime(p, r0, min_time, max_time, n_detections, t_detects_relative, number_new_symptomatic, o, rgen);
            // n_detections  = number of detected cases before (possibly early) termination 
            // t_detects = times of the detections (relative to the first detected case)
            // number_new_symptomatic = number of cases becoming symptomatic at any given timestep
            // o = outbreak structure:
            // o.individuals = list of  

            // Convert number_new_symptomatic into a list of infected individuals (day_symptomatic to day_symptomatic + infection_length-1, inclusive)
            std::vector<int> total_active_infected = number_new_symptomatic;
            MakePopulationSize(p, number_new_symptomatic, total_active_infected);
            if (p.verb == 1)
            {
                OutputPopulationDetails(p, number_new_symptomatic, total_active_infected, t_detects_relative, o);
            }

            EvaluateOutbreak(p, r0val, t_detects_relative, timepoints, total_active_infected, detections, o, results); // Is this consistent with the data?
            //std::cout << "last_time_simulated " << o.last_time_simulated << " n_t_detects " << t_detects_relative.size() << "\n"; 

        }
    }

    //Output extra statistics
    OutputAcceptanceRates (p,timepoints,results);
    if (p.more_stats==1) {
        OutputOriginTimes (p,timepoints,results);
        OutputPopulationSizes (p,timepoints,results);
        OutputProbabilityEnded (p,timepoints,results);
    }

    // Calculate statistics
    // Probability of death
    OutputOutbreakDeathStatistics(p, timepoints, results);

    // Distribution of time of initial infection
    OutputOutbreakTimeStatistics(p, timepoints, results);

    // Distribution of population size if non-zero
    OutputOutbreakPopulationStatistics(p, timepoints, results);

    return 0;
}
