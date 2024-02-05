#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>


void ConstructResults (vector<int>& timepoints, vector< vector<output> >& results);

void MakeIndexCase (run_params& p, int& t, vector<int>& t_detects, outbreak& o, gsl_rng *rgen);
void MakeIndexCaseFaster (run_params& p, int& t, const int& N, vector<int>& t_detects, const int& pp, long long& r, const long long& r_orig, vector<int>& new_incubate, vector<int>& new_detect, outbreak& o, gsl_rng *rgen);


void MakeNewCase (run_params& p, int by, vector<int>& t_detects, outbreak& o, gsl_rng *rgen);
void MakeNewCaseFaster  (run_params& p, int by, vector<int>& t_detects, const int& N, const int& pp, long long& r, const long long& r_orig, vector<int>& new_infect, vector<int>& new_incubate, vector<int>& new_detect, outbreak& o);
void SetupOutbreak(outbreak& o);
void RunSimulationTime (run_params& p, int min_time, int max_time, int& exclude, const int& n_detections, const int& N, const int& pp, long long& r, const long long& r_orig, vector<int>& new_infect, vector<int>& new_incubate, vector<int>& new_detect, vector<int>& new_number, vector<int>& t_detects, vector<int>& pop_size, outbreak& o, gsl_rng *rgen);
void MakeRelativeTime (vector<int>& t_detects, vector<int>& t_detects_relative, outbreak& o);
int CheckTermination (run_params& p, int& t, int& zeros, int& min_time, int& max_time, const int& n_detections, vector<int>& t_detects, outbreak& o);

void ConstructSummaryData (vector<int>& t_detects, vector<detect>& sim_data);


void MakePopulationSize (run_params& p, vector<int>& pop_size, vector<int>& pop_sum);
void EvaluateOutbreak (run_params& p, int& exclude, int& r0val, vector<int>& t_detects, vector<int>& timepoints, vector<int>& pop_sum, vector<detect>& detections, outbreak& o, vector< vector<output> >& results);

void CalculateAcceptance (run_params& p, int i, const vector< vector<output> >& results, vector<double>& acceptance);
