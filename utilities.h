#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>


void ConstructResults (const std::vector<int>& timepoints, std::vector< std::vector<output> >& results);


void MakeNewCase (const run_params& p, const int by, std::vector<int>& t_detects, outbreak& o, gsl_rng *rgen);

void SetupOutbreak(outbreak& o);
void RunSimulationTime (const run_params& p, const int min_time, const int max_time, const int n_detections, std::vector<int>& t_detects, std::vector<int>& pop_size, outbreak& o, gsl_rng *rgen);
void MakeRelativeTime (const std::vector<int>& t_detects, std::vector<int>& t_detects_relative, outbreak& o);
int CheckTermination (const run_params& p, int t, int zeros, int min_time, int max_time, const int n_detections, const std::vector<int>& t_detects_relative, outbreak& o);

void ConstructSummaryData (const std::vector<int>& t_detects, std::vector<detect>& sim_data);


void EvaluateOutbreak (const run_params& p, int exclude, const int r0val, std::vector<int>& t_detects, std::vector<int>& timepoints, std::vector<int>& pop_sum, std::vector<detect>& detections, outbreak& o, std::vector< std::vector<output> >& results);
void MakePopulationSize (const run_params& p, const std::vector<int>& pop_size, std::vector<int>& pop_sum);

void CalculateAcceptance (const run_params& p, const int i, const std::vector< std::vector<output> >& results, std::vector<double>& acceptance);
