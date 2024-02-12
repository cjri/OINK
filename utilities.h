#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>


void ConstructResults (const run_params& p, const std::vector<int>& timepoints, std::vector< std::vector<output> >& results);
inline void MakeNewCase (const run_params& p, const int by, std::vector<int>& t_detects, outbreak& o, gsl_rng *rgen);
void SetupOutbreak(outbreak& o);
void RunSimulationTime (const run_params& p, double r0, const int first_evaluation_time, const int last_evaluation_time, const unsigned long n_detections, std::vector<int>& t_detects, std::vector<int>& number_new_symptomatic, outbreak& o, gsl_rng *rgen);
void MakeRelativeTime (const std::vector<int>& t_detects, std::vector<int>& t_detects_relative);
int CheckTermination (const run_params& p, int t, int zeros, const int first_evaluation_time, const int last_evaluation_time, const unsigned long n_data_detections, const std::vector<int>& t_detects_relative, const outbreak& o);
void ConstructSummaryData (const std::vector<int>& t_detects, std::vector<detect>& sim_data);
void EvaluateOutbreak (const run_params& p, const unsigned long r0val, std::vector<int>& t_detects, std::vector<int>& timepoints, std::vector<int>& total_active_infected, std::vector<detect>& detections, outbreak& o, std::vector< std::vector<output> >& results);
void MakePopulationSize (const run_params& p, const std::vector<int>& number_new_symptomatic, std::vector<int>& total_active_infected);
void CalculateAcceptance (const run_params& p, const int i, const std::vector< std::vector<output> >& results, std::vector<double>& acceptance);
unsigned long CompareWithData( const std::vector<int> &timepoints, const std::vector<detect> &sim_data, const std::vector<detect> &detections);
