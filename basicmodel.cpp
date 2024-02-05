#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>


#include "basicmodel.h"
#include "io.h"
#include "pseudorandom.h"
#include "utilities.h"


int main(int argc, const char **argv) {
    //Parallels
//    MPI_Init(&argc, &argv);
//    int rank, size;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
	int seed=(int) time(NULL);
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);
	
	run_params p;
	GetOptions(p,argc,argv);
    /*ofstream case_file;
    case_file.open("Case_trajectories.dat");*/
    
    //Import data.  Currently number of cases per day.  Could do per week or coarser resolution?
    //File Detections.dat
    int n_detections=0;
    vector<detect> detections;
    ImportDetections(n_detections,detections);
    
    //We also want to know what time it is, or at least a vector of times at which to model what we know.
    int min_time;
    int max_time;
    vector<int> timepoints;
    ImportTimePoints(min_time,max_time,timepoints);
        
    //We might also want to import complex information about how the probability of detection changes each day or week.  To add here.
    
    //Next step - look at the outbreak generator.  Can stop when we have seen too many cases: This might be sooner than previous criterion.
    
    //Set up for pseudorandom code
    vector< vector<int> > all_roots;
    int N=1000000;
    vector<int> prs;
    //Vectors contain distribution outcomes
    vector<int> new_infect(N);
    vector<int> new_incubate(N);
    vector<int> new_detect(N,0);
    vector<int> new_number(N,0);
    if (p.run_fast==1) { //Use 10 different prime numbers and associated roots to make pseudorandom permutations
        GetRoots (N,prs,all_roots);
        SetupRandomParameters(p,new_infect,new_incubate,new_detect);
    }
    
    //Permutation parameters
    int pp;
    long long r;
    long long r_orig;
    if (p.run_fast==1) {
        NewPermutation (pp,r,r_orig,prs,all_roots,rgen);
    }

    
    //Parallel - could add functionality here...
//    const int num_reps = p.replicas;  // Number of replicas NB would need to make new value here
//    const int values_per_process = num_values/size;  // Divide work among processes
    // Start and end indices for each process
//    int start_index = rank * values_per_process;
//    int end_index = start_index + values_per_process;
   
    
    //Generate array of output values by R0 from 0.1 to 4.0
    //For each of these consider each time point.
    vector< vector<output> > results;
    ConstructResults(timepoints,results);
    for (p.r0=0.1;p.r0<=4.01;p.r0=p.r0+0.1) {
        if (p.run_fast==1) {
            SetupRandomPoissonParameter (p,new_number);
        }
        int r0val=floor((p.r0+0.001)*10);
        cout << "Generating simulations R0= " << p.r0 << " " << "\n";
        for (int calc=0;calc<p.replicas;calc++) { //Generate a number of simulated outbreaks at given R0
            //Construct outbreak shape
            outbreak o;
            SetupOutbreak(o);
            //Set up permutation parameters
            if (p.run_fast==1) {
                NewPermutation (pp,r,r_orig,prs,all_roots,rgen);
            }
            int exclude=0; //Flag to exclude a run
            vector<int> t_detects; //Detection times
            vector<int> pop_size;
            RunSimulationTime(p,min_time,max_time,exclude,n_detections,N,pp,r,r_orig,new_infect,new_incubate,new_detect,new_number,t_detects,pop_size,o,rgen);
            //Convert population size
            vector<int> pop_sum=pop_size;
            MakePopulationSize (p,pop_size,pop_sum);
            if (p.verb==1||p.test==1) {
                OutputPopulationDetails (p,pop_size,pop_sum,t_detects,o);
            }
            EvaluateOutbreak (p,exclude,r0val,t_detects,timepoints,pop_sum,detections,o,results);  //Is this consistent with the data?
        }
    }

    //Output extra statistics
    if (p.more_stats==1) {
        OutputAcceptanceRates (timepoints,results);
        OutputOriginTimes (timepoints,results);
        OutputPopulationSizes (timepoints,results);
        OutputProbabilityEnded (timepoints,results);
    }
    
    //Calculate statistics
    //Probability of death
    OutputOutbreakDeathStatistics (p,timepoints,results);
        
    //Distribution of time of initial infection
    OutputOutbreakTimeStatistics (p,timepoints,results);
    
    //Distribution of population size if non-zero
    OutputOutbreakPopulationStatistics (p,timepoints,results);
    
    return 0;

    

}

