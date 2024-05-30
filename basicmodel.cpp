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
//#include <omp.h>


int main(int argc, const char **argv)
{

    static gsl_rng *rgen;
    #pragma omp threadprivate(rgen)
  
    run_params p;
    GetOptions(p, argc, argv); // Fills run_params structure p with the options for this run

    // Initialize RNG to a fixed seed
    int seed = p.seed;
    gsl_rng_env_setup();
    
    //  #pragma omp parallel
    {
    rgen = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rgen, seed); // + omp_get_thread_num());
    }
    
    // Import detected infection file
    int n_detections = 0;
    std::vector<detect> detections;
    // Read Data/detections.dat
    ImportDetections(p, n_detections, detections); 
    // detections is a list of detections, with day =integer day index, cases=number of cases on that day.
    // n_detections is the sum of the number of cases.
    
    // Import observation timepoints
    int min_time;
    int max_time;
    std::vector<int> timepoints;
    // Read Data/Time_points.dat
    ImportTimePoints(p, min_time, max_time, timepoints);
    // min_time = minimum of timepoints
    // max_time = maximum of timepoints 
    // timepoints = list of integer timepoints (ordered as in file)

    // Generate array of output values by R0 from 0.1 to p.max_R0
    // For each of these consider each time point.
    std::vector<std::vector<output> > results;
    ConstructResults(p, timepoints, results);
    
    // Find extreme infection time: How long before we can conclude an outbreak has died out?
    int extreme_infection_time=floor(gsl_cdf_weibull_Pinv(0.99999,p.infection_b, p.infection_a)+1);
    
    // For each timepoint, make a list of p.R0_vals.size() zero-initialized output objects
    // Each thread has a local copy of rgen, and writes only to it's own
    // entry in the common vector results. Other variables written to are local to
    // the loop
    //#pragma omp parallel for 
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
            std::vector<int> t_detects_relative; // Detection times relative to first detected case
            std::vector<int> number_new_symptomatic; // Number of individuals becoming newly symptomatic on each day (relative to index case)

            RunSimulationTime(p, r0, min_time, max_time, n_detections, extreme_infection_time, t_detects_relative, number_new_symptomatic, o, rgen);
            // n_detections  = number of detected cases before (possibly early) termination
            // t_detects = times of the detections (relative to the first detected case)
            // number_new_symptomatic = number of cases becoming symptomatic at any given timestep
            // o = outbreak structure:

            // Convert number_new_symptomatic into a list of infected individuals (day_symptomatic to day_symptomatic + infection_length-1, inclusive)
            std::vector<int> total_active_infected = number_new_symptomatic;
            MakePopulationSize(p, o, total_active_infected);
            if (p.verb == 1)
            {
	      // Note that all this output will be interleaved in parallel version
                OutputPopulationDetails(p, number_new_symptomatic, total_active_infected, t_detects_relative, o);
            }

            EvaluateOutbreak(p, r0val, t_detects_relative, timepoints, total_active_infected, detections, o, results); // Is this consistent with the data?

        }
    } // End of parallel block

    //Output extra statistics
    OutputAcceptanceRates (p,timepoints,results);
    if (p.more_stats==1) {
        OutputOriginTimes (p,timepoints,results);
        OutputPopulationSizes (p,timepoints,results);
        OutputProbabilityEnded (p,timepoints,results);
    }

    // Calculate statistics
    // Probability that the outbreak has died out on the observation day
    OutputOutbreakDeathStatistics(p, timepoints, results);

    // Distribution of time of initial infection
    OutputOutbreakTimeStatistics(p, timepoints, results);

    // Distribution of number of actively infected individuals on the observation day
    OutputOutbreakPopulationStatistics(p, timepoints, results);

    return 0;
}
