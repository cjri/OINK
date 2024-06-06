#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

/**
 * Parses command line arguments to configure simulation parameters.
 *
 * This function iterates through the command line arguments provided to the
 * program, setting simulation parameters based on the options encountered.
 *
 * @param p A reference to a `run_params` structure where the simulation
 * parameters are stored.
 * @param argc The number of command line arguments, including the program name.
 * @param argv An array containing the command line arguments.
 *
 * Supported command line options:
 * - `--max_infections` Sets the maximum number of infections.
 * - `--replicas` Sets the number of replicas of the simulation (per R0 value).
 * - `--detect` Sets the probability of detecting an infection.
 * - `--first_detect` Sets the probability of detecting the first (index) case (by default zero)
 * - `--verb` Sets the verbosity level of the output.
 * - `--max_R0` Sets the maximum basic reproduction number (R0).
 * - `--more_stats` Enables or disables the generation of additional statistics.
 * - `--no_limit` Removes the limit on the population size growth at each timestep
 * - `--output_prefix` Sets the prefix for output files.
 * - `--input_prefix` Sets the prefix for input files.
 */
void GetOptions (run_params& p, int argc, const char **argv);

/**
 * Imports detection data from a file and populates a vector of detection records.
 *
 * This function reads detection data from a file named "Detections.dat", prefixed
 * by the `input_prefix` specified in the `run_params` structure. Each line of the
 * file is expected to contain two integers: the day of the detection and the number
 * of cases detected on that day.
 *
 * @param p A `run_params` structure containing the input file prefix.
 * @param n_detections An integer that will be updated with the total
 * number of cases detected.
 * @param detections A  vector of `detect` structures. Each `detect`
 * structure contains a day number and the number of cases detected on that day. This
 * vector will be populated with the lines read from the file.
 */
void ImportDetections(const run_params& p, int& n_detections, std::vector<detect>& detections);

/**
 * Imports time point data from a file, validates, sorts, and updates minimum and maximum time points.
 *
 * This function reads time point data from a file named "Time_points.dat", prefixed
 * by the `input_prefix` specified in the `run_params` structure.
 *
 * @param p A `run_params` structure containing the input file prefix.
 * @param min_time An integer that will be updated with the minimum
 * time point found in the file.
 * @param max_time An integer that will be updated with the maximum
 * time point found in the file.
 * @param timepoints A vector of integers that will be populated with
 * the time points read from the file.
 */
void ImportTimePoints(const run_params& p, int& min_time, int& max_time, std::vector<int>& timepoints);

/**
 * Outputs summary data of simulations to standard output.
 *
 * @param sim_data A vector of `detect` structures containing
 * the simulation data to be output.
 */
void OutputSummaryData (const std::vector<detect>& sim_data);

/**
 * Outputs detailed population statistics and detection times.
 *
 * This function prints various statistics including the numbers of new symptomatic cases on each day,
 * the total number of active infections, and details about the times of the detections relative to the first detected.
 *
 * @param p A  `run_params` structure containing simulation parameters.
 * @param number_new_symptomatic A vector containing the daily number of new symptomatic cases.
 * @param total_active_infected A vector containing the daily total number of active infections.
 * @param t_detects_relative A vector containing the timings of detections relative to the first detection.
 * @param o A `outbreak` structure containing overall outbreak details.
 *
 * The output includes:
 * - The total count and daily numbers of new symptomatic cases.
 * - The total count and daily numbers of active infections.
 * - The total number of cases and detections, along with the first detection timing and its relation to the onset of symptoms.
 * - The timings of detections relative to the first detection, indicating the spread over time.
 */
void OutputPopulationDetails (const run_params& p, const std::vector<int>& number_new_symptomatic, const std::vector<int>& total_active_infected, const std::vector<int> t_detects, outbreak& o);

/**
 * Outputs raw simulation results for specified R0 values at given timepoints.
 *
 * @param r0val The R0 value index used to select the appropriate results from the vector.
 * @param timepoints A vector of osberved timepoints
 * @param results A nested vector of `output` structures of size (timepoints.size(), p.R0_vals.size())
 */
void OutputRawData (int& r0val, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);

/**
 * Outputs the acceptance rates at different observation timepoints for multiple R0 values to files.
 *
 * This function generates a series of files, p.output_prefix + "Acceptance_rate" + (R0_index+1) + ".dat", each corresponding to a different R0 value from
 * the run parameters. Each file contains the acceptance rate P(accepted at t0 | R0), for a single R0, at specified timepoints.
 *
 * @param p The run parameters, including the output file prefix.
 * @param timepoints The timepoints at which acceptance rates are calculated.
 * @param results A nested vector of `output` structures of size (timepoints.size(), p.R0_vals.size())
 */
void OutputAcceptanceRates (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);

/**
 * Outputs the outbreak origin times for accepted simulations at different observation timepoints for multiple R0 values to files.
 *
 * This function creates a file for each R0 value specified in the run parameters, 
 * with filenames p.output_prefix + "Origin_times" + (R0_index+1)+ ".dat"
 * containing the outbreak origin times for all accepted simulations at each observation timepoint.
 * Output format consists of the timepoint followed by
 * all origin times associated with that timepoint and R0 value.
 *
 * @param p The run parameters, including the list of R0 values and the output file prefix.
 * @param timepoints The timepoints at which origin times are recorded.
 * @param results A nested vector of `output` structures of size (timepoints.size(), p.R0_vals.size())
 */
void OutputOriginTimes (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);

/**
 * Outputs the number of active infections of accepted simulations at different observation timepoints for multiple R0 values to files.
 *
 * This function generates a set of files, one for each R0 value defined in the run parameters,
 * with filenames  p.output_prefix + "Current_size" + (R0_index+1) + ".dat",
 * detailing the number of active infections of accepted simulations at various observation timepoints,
 * at the observation timepoints themselves.
 * Each file corresponds to a different R0 value; each line is the observation timepoint, and a list of all the
 * non-zero active infection counts at that observation timepoint in the accepted simulations.
 *
 * @param p The run parameters, which include the list of R0 values and the prefix for the output files.
 * @param timepoints A vector of integers indicating the timepoints at which population sizes are to be recorded.
 * @param results A nested vector of `output` structures of size (timepoints.size(), p.R0_vals.size())
 */
void OutputPopulationSizes (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);

/**
 * Outputs the proportion of accepted simulations in which the outbreak has died off before the observation time, for each R_0 value.
 *
 * This function writes to a set of files, with one file for each R0 value specified in the run parameters and
 * filenames p.output_prefix + "End_prob" + (R0_index+1) + ".dat",
 * the probability that, for an accepted simulation, the outbreak has ended prior to the observation time.
 * This is P(ended before t_o|R_0, accepted at t_o)
 *
 * @param p The run parameters, including the list of R0 values and the prefix for the output files.
 * @param timepoints A vector of integers representing the timepoints at which the probabilities are calculated.
 * @param results A 2D vector of `output` structures (timepoint_idx, R0_idx).
 */
void OutputProbabilityEnded (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);

/**
 * Outputs the probability of an outbreak died off before the observation time, averaged over all R_0 values.
 *
 * This function writes to a file, with filename p.output_prefix + "P_Outbreak_End.dat",
 * for each timepoint t_o
 * the probability that, for an accepted simulation, the outbreak has ended prior to the observation time.
 * This is P(ended before t_o|accepted at t_o)
 *
 * @param p The run parameters, which include the list of R0 values and the prefix for the output file.
 * @param timepoints A vector of integers indicating the timepoints at which the statistics are recorded.
 * @param results A 2D vector of `output` structures (timepoint_idx, R0_idx).
 */
void OutputOutbreakDeathStatistics (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);

/**
 * Outputs statistics on the date of outbreak origin at different acceptance timepoints.
 *
 * This function calculates and outputs the distribution of outbreak origin times for all simulations accepted
 * at each specified acceptance timepoint, across all R0 values provided in the run parameters.
 * Output to files with filenames p.output_prefix + "Origin_stats_" + timepoint (as a number) + ".dat"
 *
 * @param p The run parameters, including the list of R0 values and the prefix for output files.
 * @param timepoints The timepoints at which outbreak time statistics are calculated.
 * @param results A 2D vector of `output` structures (timepoint_idx, R0_idx).
 */
void OutputOutbreakTimeStatistics (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);

/**
 * Outputs the distribution of the number of active infections  for an outbreak simulation at specified time points.
 * For each time point to, this function calculates the distribution of population sizes at to across all R0 values 
 * for those simulations accepted at the observation times to and writes the distributions to separate data files
 * with filenames p.output_prefix + "Sizes_time_" + timepoint(as a number) + ".dat"
 *
 * @param p The run parameters, including the output prefix and the range of R0 values to consider.
 * @param timepoints The time points at which to output the population statistics.
 * @param results A 2D vector of `output` structures (timepoint_idx, R0_idx).
 * */
void OutputOutbreakPopulationStatistics (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results);




