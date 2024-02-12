#include "basicmodel.h"
#include "io.h"
#include "utilities.h"
#include <algorithm>
#include <string>

// Currently code assumes first infection datapoint is at time 0.

void ConstructResults(const run_params &p, const std::vector<int> &timepoints, std::vector<std::vector<output>> &results)
{
/** 
 * Initializes the results data structure for storing the outcomes of simulations across different R0 values
 * and time points.
 * 
 * Creates a vector of vectors, of size [timepoints.size, p.R0_vals.size()], with each `output` object initialized to zero.
 * 
 * @param p A constant reference to a `run_params` struct containing simulation parameters.
 * @param timepoints A constant reference to a vector of integers specifying the time points at which the
 *        outbreak is evaluated.
 * @param results A reference to a 2D vector of `output` objects where the simulation results will be stored.
 *        This vector is cleared if non-empty, and for each timepoint
 *        a vector of p.R0_vals.size() `output` objects is generated.
 */
    results.clear();
    for (unsigned long i = 0; i < timepoints.size(); i++)
    {
        std::vector<output> out;
        for (unsigned long j = 0; j < p.R0_vals.size(); j++)
        {
            output oo;
            oo.tested = 0;
            oo.accepted = 0;
            oo.dead = 0;
            out.push_back(oo);
        }
        results.push_back(out);
    }
}

inline void MakeNewCase(const run_params &p, const int by, std::vector<int> &t_detects, outbreak &o, gsl_rng *rgen)
{
    /**
     * @brief Creates a new case of infection within a given outbreak simulation.
     *
     * This function generates a new patient case. It calculates the time of infection and symptom onset based on Weibull distributions, with parameters stored in `p`, determines randomly whether the case is detected, and updates the outbreak data structure accordingly.
     * The detection times of new cases (if detected) are added to the `t_detects` vector, and the outbreak's `time_first_detect` is updated if this case is detected earlier than any previously detected case.
     * 
     * @param p A reference to a `run_params` structure containing parameters for the infection and incubation periods, detection probabilities, and the time from symptom onset to detection.
     * @param by A value of -1 indicates the case is a primary (index) case not infected by anyone within the simulation.
     * @param t_detects A reference to a vector of integers storing the times (relative to the start of the simulation) at which cases were detected. If the new cases is detected the detection time is appended to this list
     * @param o A reference to the outbreak structure where the new case will be added. 
     * @param rgen A pointer to a GSL random number generator.
     * 
     */

    patient pt;
    // Time of infection
    if (by == -1) // Index case
    {
        pt.time_infected = 0;
    }
    else
    {
      pt.time_infected = o.individuals[by].time_symptom_onset + static_cast<int>(floor(gsl_ran_weibull(rgen, p.infection_b, p.infection_a)+0.5)); 
    }
    pt.time_symptom_onset = pt.time_infected + static_cast<int>(floor(gsl_ran_weibull(rgen, p.incubation_b, p.incubation_a)+0.5));
    if (by == -1) // Index case
    {
        pt.detected = gsl_ran_bernoulli(rgen, p.probability_first_detect);
    }
    else
    {
        pt.detected = gsl_ran_bernoulli(rgen, p.probability_detect);
    }
    if (pt.detected == 1)
    {
        pt.time_reported = pt.time_symptom_onset + p.time_symptom_onset_to_detect; 
        if (o.time_first_detect == -1 || pt.time_reported < o.time_first_detect)
        {
            o.time_first_detect = pt.time_reported;
        }
        t_detects.push_back(pt.time_reported);
    }
    else
    {
        pt.time_reported = -1;
    }
    o.individuals.push_back(pt);
}

void SetupOutbreak(outbreak &o)
{
    /* Initialize an outbreak structure */
    o.time_first_detect = -1; 
    o.last_time_simulated = 0;
}

void RunSimulationTime(const run_params &p, 
                       const double r0,
                       const int first_evaluation_time, 
                       const int last_evaluation_time, 
                       const unsigned long n_data_detections, 
                       std::vector<int> &t_detects_relative, 
                       std::vector<int> &number_new_symptomatic, // Population size at integer time t, relative to start of the simulation
                       outbreak &o, 
                       gsl_rng *rgen)
{
/**
 * Main simulation function
 * 
 * This function simulates the spread of an outbreak from an index case over a predefined period, or until certain
 * termination conditions are met. It tracks the number of new infections and symptomatic individuals at each time
 * step, applying a limit to the number of new cases generated each day to maintain computational efficiency.
 * 
 * @param p A `run_params` struct containing simulation parameters.
 * @param r0 The reproduction number
 * @param first_evaluation_time The earliest time to evaluate whether the simulations match the data (timepoints). 
 * @param last_evaluation_time The last timepoint
 * @param n_detections The total number of detections in the data.
 * @param t_detects Detection times (relative) to the first detected case ; during the function these are absolute detection times
 * @param number_new_symptomatic The number of new infections at each integer time t, time measured relative to the start of the simulation. 
 * @param o An `outbreak` object
 * @param rgen A pointer to a GSL random number generator
 * 
 * @note Detection times are adjusted to be relative to the first detection, and the simulation may terminate early if
 *   it meets certain criteria, such as a prolonged period without new cases or reaching a maximum number of detections.
 */

    int t = 0;
    std::vector<int> t_detects;
    MakeNewCase(p, -1, t_detects, o, rgen); // Make index case
    int index = -1; // 
    int zeros = 0; // Count of the number of timesteps with no new cases
    while (t < p.max_simulation_time)
    {
        index++;
        int added = 0;
        number_new_symptomatic.push_back(0); // Count of the number of people infected at exactly time t
        for (unsigned long i = 0; i < o.individuals.size(); i++)
        { 
            if (int(o.individuals[i].time_infected) == t) 
            {   
                number_new_symptomatic[index]++; 
            }
            if (int(o.individuals[i].time_symptom_onset) == t) 
            { // Generate infections for cases which become symptomatic on this day
                zeros = 0; // Reset zeros 
                int infect = 0; // How many new cases generated by this individual
                if (added < p.add_limit)
                { // We have a limit here on how many cases to generate each day.
                    // The limit can be effectively removed with the --no_limit flag.
                    // Having a limit makes the code run faster for outbreaks with larger R0.  The default limit of 1000 ensures that some cases
                    // will almost certainly be detected each day, therefore altering the stats.
                    // Some thought needs to go into what happens if there are lots of detections each day in the data: Probably that isn't a good application of this code.
                    infect = gsl_ran_poisson(rgen, r0);
                    added += infect;
                }
                for (int j=0;j<infect;j++) {  
                    //Add infect new cases of infection
                    MakeNewCase (p,i,t_detects,o,rgen); /* These should all be infected at new time >=t */
                }
            }
        }
        zeros++;
        o.last_time_simulated = t;

        // Check whether to terminate the simulation (early)
        std::sort(t_detects.begin(), t_detects.end());
        int term = CheckTermination(p, t, zeros, first_evaluation_time, last_evaluation_time, n_data_detections, t_detects, o);

        if (term == 1)
        {
            break;
        }
        /* What's consistent at this point (time t):
        we've counted all the individuals infected <= t
        added new infections for those individuals who are *also* symptomatic <=t
        t_detect is complete up to t, but contains a subset of 
        detections at times >t.
        zeros is 1 if a case became symptomatic in this timestep;
        so is the number of days *plus 1* since the last symptomatic case.
        Have possibly decided to terminate early.
        */

        t++;
    }

    // Calculate relative detection times
    MakeRelativeTime(t_detects, t_detects_relative);
    // Ensure relative detection times are sorted (Should be from l 139)
    std::sort(t_detects_relative.begin(), t_detects_relative.end()); 
}

void MakeRelativeTime(const std::vector<int> &t_detects, 
                      std::vector<int> &t_detects_relative)
{
/**
 * Transforms absolute detection times into relative times with respect to the earliest detection.
 * 
 * This function takes a vector of detection times relative to the start of the simulation (`t_detects`), and ouputs a vector
 * `t_detects_relative` of these times relative to the first detection.
 * 
 * @param t_detects The times of detections.
 * @param t_detects_relative Relative detection times,
 */

    t_detects_relative.clear();
    if(t_detects.empty()) 
    {
        return;
    }

    int min = std::numeric_limits<int>::max();
    // Find min(t_detects)
    for (unsigned long i = 0; i < t_detects.size(); i++)
    {
        if (t_detects[i] < min)
        {
            min = t_detects[i];
        }
    }
    // If sorted could use min=t_detects[0]

    t_detects_relative.clear();
    t_detects_relative.reserve(t_detects.size());
    // Subtract min fom all detection times
    for (unsigned long i = 0; i < t_detects.size(); i++)
    {
        t_detects_relative.push_back(t_detects[i] - min);
    }
    // Fix the time of the index infection
}

int CheckTermination(const run_params &p, 
                     const int t, 
                     const int zeros, 
                     const int first_evaluation_time, 
                     const int last_evaluation_time,
                     const unsigned long n_detections, 
                     const std::vector<int> &t_detects, 
                     const outbreak &o)
{
/**
 * Determines whether the simulation of an outbreak should be terminated based on specific criteria.
 * 
 * This function checks for several conditions to decide if the simulation of an outbreak should be terminated.
 * It considers the number of days with zero new symptoms, the total number of detections relative to the dataset, and
 * the timing of these detections against predefined parameters. Termination conditions include: no new symptoms for a
 * period longer than the infection length plus one, and some matching conditions TODO.
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



    // Terminate simulation if there have been more than infection_length days with no individuals developing symptoms on that day (zeros is the number of days +1 since an individual last developed symptoms)
    // With Flu parameters this is probably OK, but it will require time to symptoms and time to infection to be almost always less than infection_length, which may not be the case for other parameter choices.
    if (zeros > p.infection_length + 1)
    {
        if (p.verb == 1)
        {
            std::cout << "Terminate: Outbreak has died out\n";
        }
        // cout << "Zero out\n";
        return 1;
    }

    if ((o.time_first_detect != -1) && (t - o.time_first_detect >= last_evaluation_time)) 
    {
        // This stops the simulation once the time has reached the last observation. 
        // (At this point t >= time_first_detect, so first detection is known)
        // Matching will not consider any cases after this time
        if (p.verb == 1)
        {
            std::cout << "Terminate: End of time after last detection\n";
        }
        return 1;
    }  



    if (t_detects.size() > n_detections)
    {
        if ((o.time_first_detect != -1) && (t  >= t_detects[n_detections])) 
        {
            // At this point we know that all timepoints on or after t_detects[n_detections]
            // will fail, and all those before are completely simulated.
            if(p.verb==1)
            {
                std::cout << "Terminate: End of time after last detection\n";
            }
            return 1;
        }
    }

    // If we have simulated up to the first detection, are there too many cases before the first observation timepoint.
    // In this case, we know that all the matching comparisons with the data will fail.
    // Need to have simulated up to first detection, otherwise simulation could generate a new detected case in (t, t_detect[0]]
    // which might be valid for [t_new_detect, t_new_detect+first_evaluation_time]
    if((o.time_first_detect != -1) && (t >= o.time_first_detect)) {
        // Test whether the number of detections on or before the first evaluation time-point is greater than the number of detections in the dataset - if so, this will be the case for all later ones.
        unsigned long det_early = 0;
        // Now assume sorted - break if time greater than detection time
        for (unsigned long i = 0; (i < t_detects.size()) && (t_detects[i]<=o.time_first_detect + first_evaluation_time); i++)
        {
            // Number of detections happening on or before the first evaluation time-point
            det_early++;
        }

        if (det_early > n_detections)
        { // This number exceeds the total number of detections: Easy flag for non-compatibility with the data
            if (p.verb == 1)
            {
                std::cout << "Terminate: Too many detections too early\n";
            }
            return 1;
        }
    }
    return 0;
    
}


void ConstructSummaryData(const std::vector<int> &t_detects, std::vector<detect> &sim_data)
{

/**
 * Converts a list of detection times from a simulation into the same format as the input data.
 * 
 * This function takes a sorted vector of detection times (`t_detects`) and appends it a vector of `detect` objects (`sim_data`),
 * where each `detect` object contains the day relative to the first detected case and the number of cases detected on that day. 
 * Most of the time, it is likely that you want to pass an empty vector to this function.
 * If t_detects is empty, it creates a single detection object with day=0 and cases=0, which is important for the
 * current implementation of the matching function.
 * If t_detects.size()>0 and t_detects[0] !=0 then it also inserts this detection object; however, for all inputs currently
 * t_detects[0]=0
 * 
 * @param t_detects A sorted vector of integers representing the times (days) at which detections occurred during the simulation,
 *                  relative to the first detected case. 
 * @param sim_data A vector of `detect` objects to which the summary of detections will be appended.
 */

    if(t_detects.empty()) {
        return;
    }
    // Produces counts describing detections in the simulation in a format that matches the input data
    // Assume t_detects is already sorted
    int index = 0;
    int count = 0;
    //std::sort(t_detects.begin(), t_detects.end());
    for (unsigned long i = 0; i < t_detects.size(); i++)
    {
        if (t_detects[i] == index)
        {
            count++;
        }
        else
        {
            detect d;
            d.day = index;
            d.cases = count;
            sim_data.push_back(d);
            index = t_detects[i];
            count = 1;
        }
    }
    detect d;
    d.day = index;
    d.cases = count;
    sim_data.push_back(d);
}

unsigned long CompareWithData( const std::vector<int> &timepoints, const std::vector<detect> &sim_data, const std::vector<detect> &detections)
{
    // Note that if timepoint[i] does not match, all later timepoints also do not match.
    // So we can count the number of accepted timepoints
    unsigned long n_timepoints_accepted=0;
    unsigned long idx_sim=0;

    while (n_timepoints_accepted<timepoints.size()) 
    { 
        int accept=1;
        while ((idx_sim<sim_data.size()) && (sim_data[idx_sim].day <= timepoints[n_timepoints_accepted]))
        {
            if ((idx_sim>=detections.size()) || (sim_data[idx_sim].day != detections[idx_sim].day) || (sim_data[idx_sim].cases != detections[idx_sim].cases))
            {
                // Mismatch with data
                accept = 0;
                break;
            }
            idx_sim++;
        }
        if((idx_sim<detections.size()) && (detections[idx_sim].day <=timepoints[n_timepoints_accepted])) {
            accept=0;
        }
        if(accept==0) 
        {
            break;
        }
        n_timepoints_accepted++;

    }
    if((n_timepoints_accepted>0) && (sim_data.size()==0))
    {
        std::cout << "ERR " << n_timepoints_accepted << " " << sim_data.size() << " " << detections.size() << "\n";
        for(auto t : detections)
            std::cout << t.day << " " << t.cases << "\n";
        for(auto t : timepoints)
            std::cout << "tp: "<< t << "\n";
        std::cout << "\n";
    }

    return n_timepoints_accepted;
}

void EvaluateOutbreak(const run_params &p, const unsigned long r0val, std::vector<int> &t_detects_relative, std::vector<int> &timepoints, std::vector<int> &total_active_infected, std::vector<detect> &detections, outbreak &o, std::vector<std::vector<output>> &results)
{
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

    for (unsigned long i = 0; i < timepoints.size(); i++)
    {                               // Cycle through time points at which we are evaluating the simulation
        results[i][r0val].tested++; // Count of number of simulations
    }
    // Check the simulation against the detection data
    std::vector<detect> sim_data;
    ConstructSummaryData(t_detects_relative, sim_data); // Could generate these from absolute times + o.first_time_detected
    if (p.verb == 1)
    {
        OutputSummaryData(sim_data);
    }

    unsigned long n_timepoints_accepted = CompareWithData(timepoints, sim_data, detections);

    for (unsigned long i=0;i<n_timepoints_accepted;i++) 
    { 
    // Simulation fits the data
        if (p.verb == 1)
        {
            std::cout << "Accepted at time point " << timepoints[i] << "\n";
            std::cout << "Size of total_active_infected " << total_active_infected.size() << "\n";
            int time = o.time_first_detect + timepoints[i];
            std::cout << "Time " << time << "\n";
            if (total_active_infected.size() > time)
            { // Current number of infections
                std::cout << "Size " << total_active_infected[time] << "\n";
            }
        }
        results[i][r0val].accepted++;                           // Record acceptance
        results[i][r0val].origin_time.push_back(-o.time_first_detect); // Origin time of outbreak relative to first detection
        // What do we do if first outbreak time is not 0 in data???
        if(o.time_first_detect<0) {
            std::cout << "Negative first detection time " << o.time_first_detect << " " << n_timepoints_accepted << " " << sim_data.size() << "\n";
        }
        // Find current size
        int time = o.time_first_detect + timepoints[i];
        if (total_active_infected.size() > time)
        { // Store current number of "active" infections at this time
            results[i][r0val].current_size.push_back(total_active_infected[time]);
        }
        else
        {
            results[i][r0val].current_size.push_back(0);
            results[i][r0val].dead++;
        }
}
    if (p.verb == 1)
    {
        std::cout << "\n";
    }
}

void MakePopulationSize(const run_params &p, const std::vector<int> &number_new_symptomatic, std::vector<int> &total_active_infected)
{
/**
 * Calculates the population size affected by infections over time.
 * 
 * This function processes an input array representing the population size at each time point (number_new_symptomatic) and spreads
 * the effect of each infection over a period defined by `p.infection_length`. For each infected individual at time
 * point `i`, their impact is added to the total population size (`total_active_infected`) from time `i` to `i + p.infection_length - 1`,
 * inclusive. 
 * 
 * @param p A `run_params` struct containing simulation parameters, specifically the length
 *          of time an infection affects the population (`infection_length`).
 * @param number_new_symptomatic The number of newly symptomatic patients at each
 *          time point.
 * @param total_active_infected Output of the current population size affected by infections

 *  @note The function assumes total_active_infected is of the correct size.
 **/

    for (unsigned long i = 0; i < total_active_infected.size(); i++)
    {
        total_active_infected[i] = 0;
    }
    for (unsigned long i = 0; i < total_active_infected.size() && i < number_new_symptomatic.size(); i++)
    {
        if (number_new_symptomatic[i] > 0)
        {
            for (int j = 0; j < p.infection_length; j++)
            {
                if (i + j < total_active_infected.size())
                {
                    total_active_infected[i + j] = total_active_infected[i + j] + number_new_symptomatic[i];
                }
            }
        }
    }
}

void CalculateAcceptance(const run_params &p, const int i, const std::vector<std::vector<output>> &results, std::vector<double> &acceptance)
{
/**
 * Calculates and normalizes the acceptance rates for different R0 values at a specific time point.
 * 
 * This function computes the acceptance rate for simulations at a given time point `i` across a range of R0 values.
 * The acceptance rate for each R0 value is calculated as the ratio of accepted simulations to the total number of
 * tested simulations. These raw acceptance rates are then normalized across all R0 values considered, so that the
 * sum of acceptance rates for all R0 values equals 1.
 * 
 * @param p A `run_params` struct containing simulation parameters.
 * @param i The index of the time point for which acceptance rates are being calculated.
 * @param results A 2D vector of `output` objects.
 * @param acceptance Output acceptance rates, normalized such that their sum equals 1.
 */

    double tot = 0.0;
    acceptance.clear();
    for (unsigned long r0val=0; r0val<p.R0_vals.size(); r0val++)
    {
        double acc = static_cast<double>(results[i][r0val].accepted ) / results[i][r0val].tested;
        tot = tot + acc;
        acceptance.push_back(acc);
    }
    for(unsigned long j = 0; j < acceptance.size(); j++)
    {
        acceptance[j] = acceptance[j] / tot;
    }
}
