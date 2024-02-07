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
    ImportDetections(n_detections, detections); 
    // detections is a list of detections, with day =integer day index, cases=number of cases on that day.
    // n_detections is the sum of the number of cases.

    // We also want to know what time it is, or at least a std::vector of times at which to model what we know.
    int min_time;
    int max_time;
    std::vector<int> timepoints;
    // Read Data/Time_points.dat
    ImportTimePoints(min_time, max_time, timepoints);
    // min_time = minimum of 0 and timepoints
    // max_time = maximum of 1000 and timepoints 
    // timepoints = list of integer timepoints (ordered as in file)

    // Generate array of output values by R0 from 0.1 to 4.0
    // For each of these consider each time point.
    std::vector<std::vector<output>> results;
    ConstructResults(timepoints, results);
    // For each timepoint, make a list of 51? zero-initialized output objects
    // Fix this loop over R0 Values
    for (p.r0 = 0.1; p.r0 <= 4.01; p.r0 = p.r0 + 0.1)
    {
        int r0val = floor((p.r0 + 0.001) * 10); // Some sort of index
        std::cout << "Generating simulations R0= " << p.r0 << " "
                  << "\n";
        for (int calc = 0; calc < p.replicas; calc++)
        { // Generate a number of simulated outbreaks at given R0
            // Construct outbreak shape
            outbreak o;
            SetupOutbreak(o); // Initialize outbreak

            // Set up permutation parameters
            int exclude = 0;            // Flag to exclude a run
            std::vector<int> t_detects; // Detection times ?
            std::vector<int> pop_size; // Population sizes ?

            RunSimulationTime(p, min_time, max_time, n_detections, t_detects, pop_size, o, rgen);
            // n_detections  = number of detected cases before (possibly early) termination 
            // t_detects = times of the detections (relative to the first detected case)
            // pop_size = number of cases becoming symptomatic at any given timestep
            // o = outbreak structure:
            // o.individuals = list of  

            // Convert pop_size into a list of infected individuals (day_symptomatic to day_symptomatic + infection_length-1, inclusive)
            std::vector<int> pop_sum = pop_size;
            MakePopulationSize(p, pop_size, pop_sum);
            if (p.verb == 1)
            {
                OutputPopulationDetails(p, pop_size, pop_sum, t_detects, o);
            }
            EvaluateOutbreak(p, exclude, r0val, t_detects, timepoints, pop_sum, detections, o, results); // Is this consistent with the data?
        }
    }

    // Output extra statistics
    if (p.more_stats == 1)
    {
        OutputAcceptanceRates(timepoints, results);
        OutputOriginTimes(timepoints, results);
        OutputPopulationSizes(timepoints, results);
        OutputProbabilityEnded(timepoints, results);
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
