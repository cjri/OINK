#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

void GetOptions (run_params& p, int argc, const char **argv);
void ImportDetections(const run_params& p, int& n_detections, std::vector<detect>& detections);
void ImportTimePoints(const run_params& p, int& min_time, int& max_time, std::vector<int>& timepoints);
//void OutputCaseData (int& first_detect, ofstream case_file, std::vector<pat>& pdat);


void OutputSummaryData (const std::vector<detect>& sim_data);
void OutputPopulationDetails (const run_params& p, const std::vector<int>& number_new_symptomatic, const std::vector<int>& total_active_infected, const std::vector<int> t_detects, outbreak& o);
void OutputRawData (int& r0val, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);
void OutputAcceptanceRates (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);
void OutputOriginTimes (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);
void OutputPopulationSizes (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);
void OutputProbabilityEnded (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);
void OutputOutbreakDeathStatistics (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);
void OutputOutbreakTimeStatistics (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);
void OutputOutbreakPopulationStatistics (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);




