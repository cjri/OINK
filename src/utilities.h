#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

/** 
 * Initializes the results data structure for storing the outcomes of simulations across different R0 values
 * and time points.
 * 
 * Creates a vector of vectors, of size [timepoints.size, p.R0_vals.size()], with each `output` object initialized to zero.
 * 
 * @param p A  `run_params` struct containing simulation parameters.
 * @param timepoints Avector of integers specifying the time points at which the
 *        outbreak is evaluated.
 * @param results A reference to a 2D vector of `output` objects where the simulation results will be stored.
 *        This vector is cleared if non-empty, and for each timepoint
 *        a vector of p.R0_vals.size() `output` objects is generated.
 */
void ConstructResults (const run_params& p, const std::vector<int>& timepoints, std::vector< std::vector<output> >& results);

/**
 * @brief Creates a new case of infection within a given outbreak simulation.
 *
 * This function generates a new patient case. It calculates the time of infection and symptom onset based on Weibull distributions, with parameters stored in `p`, determines randomly whether the case is detected, and updates the outbreak data structure accordingly.
 * The detection times of new cases (if detected) are added to the `t_detects` vector, and the outbreak's `time_first_detect` is updated if this case is detected earlier than any previously detected case.
 * 
 * @param p A `run_params` structure containing parameters for the infection and incubation periods, detection probabilities, and the time from symptom onset to detection.
 * @param by A value of -1 indicates the case is a primary (index) case not infected by anyone within the simulation.
 * @param t_detects A reference to a vector of integers storing the times (relative to the start of the simulation) at which cases were detected. If the new cases is detected the detection time is appended to this list
 * @param t_detects_after A reference to a vector of integers storing the times (relative to the start of the simulation) at which cases were detected (with increased surveillance)
 * @param o A reference to the outbreak structure where the new case will be added. 
 * @param rgen A pointer to a GSL random number generator.
 * 
 */
inline void MakeNewCase (const run_params& p, const int by, std::vector<int>& t_detects, std::vector<int>& t_detects_after , outbreak& o, gsl_rng *rgen);

    
/* Initialize an empty outbreak structure 
 * @param o A reference to the outbreak structure.
*/
void SetupOutbreak(outbreak& o);

/**
 * Main simulation function
 * 
 * This function simulates the spread of an outbreak from an index case over a predefined period, or until certain
 * termination conditions are met. It tracks the number of new infections and newly symptomatic individuals at each time
 * step, applying a limit to the number of new cases generated each day to maintain computational efficiency.
 * 
 * @param p A `run_params` struct containing simulation parameters.
 * @param r0 The reproduction number
 * @param first_evaluation_time The earliest time to evaluate whether the simulations match the data (timepoints). 
 * @param last_evaluation_time The last timepoint
 * @param n_detections The total number of detections in the data.
 * @param t_detects_relative Detection times (relative) to the first detected case 
 * @param number_new_symptomatic The number of new infections at each integer time t, time measured relative to the start of the simulation. 
 * @param o An `outbreak` object
 * @param rgen A pointer to a GSL random number generator
 **/
void RunSimulationTime (const run_params& p, const double r0, const int first_evaluation_time, const int last_evaluation_time, const unsigned long n_detections, const int extreme_infection_time, std::vector<int>& t_detects_relative, std::vector<int>& number_new_symptomatic, outbreak& o, gsl_rng *rgen);

/**
 * Transforms absolute detection times into relative times with respect to the earliest detection.
 * 
 * This function takes a vector of detection times relative to the start of the simulation (`t_detects`), and ouputs a vector
 * `t_detects_relative` of these times relative to the first detection.
 * 
 * @param t_detects The times of detections.
 * @param t_detects_relative Relative detection times,
 */
void MakeRelativeTime (const std::vector<int>& t_detects, std::vector<int>& t_detects_relative);

/**
 * Determines whether the simulation of an outbreak should be terminated (early) based on specific criteria.
 * 
 * This function checks for several conditions to decide if the simulation of an outbreak should be terminated.
 * It considers the number of days with zero new symptoms, the total number of detections relative to the dataset, and
 * the timing of these detections against predefined parameters. Termination conditions include: no new symptoms for a
 * period longer than the infection length plus one, and some more complex matching conditions.
 * 
 * @param p A `run_params` struct containing simulation parameters.
 * @param t The current simulation time. (Test is at end of time t).
 * @param zeros The number + 1 of consecutive days with zero new cases becoming symptomatic.
 * @param first_evaluation_time The first time at which to evaluate consistency with the data
 * @param last_evaluation_time The time of the last evaluation time
 * @param n_detections The number of detections specified by the dataset.
 * @param t_detects The times of detections relative to the start of the simulation
 * @param o An `outbreak` object
 * 
 * @return Returns 1 if any of the termination criteria are met, indicating that the simulation should be terminated.
 *         Otherwise, returns 0, indicating that the simulation should continue.
 * 
 * @note If the simulation terminates, we know that either 
 *       - the outbreak has died out, so our data is complete for all time
 *       - we have complete data up to all timepoints that will be tested
 *       - we have complete data up to time t, which is on or after the first detected case, and know that all (later) timepoints will fail to match
 */
int CheckTermination(const run_params &p, const int t, const int zeros, const int extreme_infection_time, const int first_evaluation_time, const int last_evaluation_time, const unsigned long n_detections, const std::vector<int> &t_detects, const outbreak &o);

/**
 * Converts a list of detection times from a simulation into the same format as the input data.
 * 
 * This function takes a sorted vector of detection times (`t_detects`) and appends to it a vector of `detect` objects (`sim_data`),
 * where each `detect` object contains the day relative to the first detected case and the number of cases detected on that day. 
 * It is likely that you want to pass an empty vector to this function.
 * If t_detects.size()>0 and t_detects[0] !=0 then it also inserts a detection with day=0 and cases=0; however, for all inputs currently
 * t_detects[0]=0 or t_detects is empty.
 * 
 * @param t_detects A sorted vector of integers representing the times (days) at which detections occurred during the simulation,
 *                  relative to the first detected case. 
 * @param sim_data A vector of `detect` objects to which the summary of detections will be appended.
 */
void ConstructSummaryData (const std::vector<int>& t_detects, std::vector<detect>& sim_data);

/**
 * Compares a sequence of timepoints with simulation data and actual detection data to count the number of timepoints 
 * where the simulation data matches the detection data up to the first mismatch.
 * 
 * This function iterates over a list of timepoints and compares them against corresponding simulation data (`sim_data`)
 * and actual detection data (`detections`). A timepoint is accepted if, for all simulation data up to that timepoint,
 * the day and case count match exactly between the simulation data and the detection data. The comparison stops at the 
 * first timepoint where a mismatch occurs, as all later timepoints by necessity do not match.
 * 
 * Note:
 * - If a timepoint in `timepoints` does not have a matching day in `sim_data` or `detections`, or if the case counts do not match,
 *   the function will stop checking further timepoints.
 * - The function assumes `sim_data` and `detections` are sorted by the day.

 * 
 * @param timepoints Timepoints to be compared.
 * @param sim_data A vector of `detect` structures representing simulation data
 * @param detections A vector of `detect` structures representing actual detection data, similar to `sim_data`.
 * @return The number of accepted timepoints where simulation data matches detection data up to the first detected mismatch.
 */
unsigned long CompareWithData( const std::vector<int> &timepoints, const std::vector<detect> &sim_data, const std::vector<detect> &detections);

/**
 * Evaluates the outcome of an outbreak simulation at specified time points against detection data.
 * 
 * This function cycles through specified time points. For each time point, it checks if the simulated data match the real detection data up to that
 * point. The evaluation considers the number of cases detected and the days on which these detections occurred.
 * The function updates results with the number of simulations tested and accepted for each R0 value and
 * time point, along with the origin time of the outbreak, the current size of the population affected, and the number of simulations for which the outbreak is over,  
 * at each accepted time point.
 * 
 * @param p A `run_params` struct containing simulation parameters.
 * @param r0val The R0 value index being evaluated, used to index into the results matrix.
 * @param t_detects_relative The times at which detections occurred in the simulation, relative to the first detected case
 * @param timepoints The time points at which the outbreak is evaluated. Sorted in increasing order.
 * @param total_active_infected The number of active infections at each time point.
 * @param detections The actual detection data for comparison.
 * @param o An `outbreak` object containing information about the outbreak being simulated.
 * @param results A reference to a vector of vectors of `output` objects of size [timepoints.size(), p.R0_vals.size()] 
 *        where the evaluation results are stored. Each element
 *        stores the number of simulations tested, accepted,
 *        the origin time of accepted outbreaks, the current size of the 
 *        outbreak at each accepted time point, and the number of simulations
 *        where the outbreak has died out
 */
void EvaluateOutbreak (const run_params& p, const unsigned long r0val, std::vector<int>& t_detects, std::vector<int>& timepoints, std::vector<int>& total_active_infected, std::vector<detect>& detections, outbreak& o, std::vector< std::vector<output> >& results);

/**
 * Calculates the population size affected by infections over time.
 * 
 * This function processes the outbreak data structure which contains (in its individuals field)
 * details about all infected individuals. For each infected individual (i) at 
 * their impact is added to the total population size (`total_active_infected`) from time `i.time_symptom_onset` to 
 * `i.time_symptom_onset + i.infection_length - 1`, inclusive. 
 * 
 * @param p A `run_params` struct containing simulation parameters, specifically the length
 *          of time an infection affects the population (`infection_length`).
 * @param o An `outbreak` object containing information about the outbreak being simulated.
 * @param total_active_infected Output of the current population size affected by infections

 *  @note The function assumes total_active_infected is of the correct size.
 **/
void MakePopulationSize (const run_params& p, const outbreak& o, std::vector<int>& total_active_infected);
