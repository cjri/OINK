#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

void GetOptions (run_params& p, int argc, const char **argv);
void ImportDetections(int& n_detections, vector<detect>& detections);
void ImportTimePoints(int& min_time, int& max_time, vector<int>& timepoints);
//void OutputCaseData (int& first_detect, ofstream case_file, vector<pat>& pdat);
void OutputSummaryData (vector<detect>& sim_data);
void OutputPopulationDetails (run_params& p, vector<int>& pop_size, vector<int>& pop_sum, vector<int> t_detects, outbreak& o);
void OutputRawData (int& r0val, vector<int>& timepoints, vector< vector<output> >& results);
void OutputAcceptanceRates (vector<int>& timepoints, vector< vector<output> >& results);
void OutputOriginTimes (vector<int>& timepoints, vector< vector<output> >& results);
void OutputPopulationSizes (vector<int>& timepoints, vector< vector<output> >& results);
void OutputProbabilityEnded (vector<int>& timepoints, vector< vector<output> >& results);
void OutputOutbreakDeathStatistics (run_params& p, vector<int>& timepoints, vector< vector<output> >& results);
void OutputOutbreakTimeStatistics (run_params& p, vector<int>& timepoints, vector< vector<output> >& results);
void OutputOutbreakPopulationStatistics (run_params& p, vector<int>& timepoints, vector< vector<output> >& results);




