#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

void GetOptions (run_params& p, int argc, const char **argv);
void ImportDetections(int& n_detections, std::vector<detect>& detections);
void ImportTimePoints(int& min_time, int& max_time, std::vector<int>& timepoints);
//void OutputCaseData (int& first_detect, ofstream case_file, std::vector<pat>& pdat);
void OutputSummaryData (std::vector<detect>& sim_data);
void OutputPopulationDetails (run_params& p, std::vector<int>& pop_size, std::vector<int>& pop_sum, std::vector<int> t_detects, outbreak& o);
void OutputRawData (int& r0val, std::vector<int>& timepoints, std::vector< std::vector<output> >& results);
void OutputAcceptanceRates (run_params& p, std::vector<int>& timepoints, std::vector< std::vector<output> >& results);
void OutputOriginTimes (run_params& p, std::vector<int>& timepoints, std::vector< std::vector<output> >& results);
void OutputPopulationSizes (run_params& p, std::vector<int>& timepoints, std::vector< std::vector<output> >& results);
void OutputProbabilityEnded (run_params& p, std::vector<int>& timepoints, std::vector< std::vector<output> >& results);
void OutputOutbreakDeathStatistics (run_params& p, std::vector<int>& timepoints, std::vector< std::vector<output> >& results);
void OutputOutbreakTimeStatistics (run_params& p, std::vector<int>& timepoints, std::vector< std::vector<output> >& results);
void OutputOutbreakPopulationStatistics (run_params& p, std::vector<int>& timepoints, std::vector< std::vector<output> >& results);




