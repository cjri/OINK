#include "basicmodel.h"
#include "io.h"
#include "utilities.h"
#include <algorithm>
#include <string>

void ConstructResults(const std::vector<int> &timepoints, std::vector<std::vector<output>> &results)
{
    /* Initialize results vector of vectors. Unclear why 51 vectors at each output timepoint */

    for (unsigned int i = 0; i < timepoints.size(); i++)
    {
        std::vector<output> out;
        for (int i = 0; i < 51; i++)
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

void MakeNewCase(const run_params &p, const int by, std::vector<int> &t_detects, outbreak &o, gsl_rng *rgen)
{
    /**
     * @brief Creates a new case of infection within a given outbreak simulation.
     *
     * This function generates a new patient case based on the provided parameters and random number generator. It assigns the patient's infector, calculates the time of infection and symptom onset based on Weibull distributions, determines whether the case is detected, and updates the outbreak data structure accordingly.
     *
     * @param p A reference to a `run_params` structure containing parameters for the infection and incubation periods, detection probabilities, and the time from symptom onset to detection.
     * @param by The index of the patient who infected the new case. A value of -1 indicates the case is a primary (index) case not infected by anyone within the simulation.
     * @param t_detects A reference to a vector of integers where the times of reported detection for detected cases will be appended.
     * @param o A reference to the outbreak structure where the new case will be added. This structure contains a vector of individual patient cases and other relevant outbreak information.
     * @param rgen A pointer to a GSL random number generator used for drawing random variates needed for simulating the infection and detection processes.
     *
     * The function modifies the `o` (outbreak) argument by adding a new `patient` to its list of individuals. The attributes of this new patient, such as the time of infection, time of symptom onset, detection status, and time reported (if detected), are calculated within the function. The detection times of new cases (if detected) are added to the `t_detects` vector, and the outbreak's `time_first_detect` is updated if this case is detected earlier than any previously detected case.
     *
     * Note: The function uses Weibull distributions for modeling the time until infection and the incubation period, and Bernoulli distributions for modeling the probability of detection. The parameters for these distributions are taken from the `run_params` structure.
     */

    patient pt;
    pt.infected_by = by;
    // Time of infection
    if (by == -1) // Index case
    {
        pt.time_infected = 0.0;
    }
    else
    {
        pt.time_infected = o.individuals[by].time_symptom_onset + gsl_ran_weibull(rgen, p.infection_b, p.infection_a); // This is time as a double
    }
    pt.time_symptom_onset = pt.time_infected + gsl_ran_weibull(rgen, p.incubation_b, p.incubation_a);
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
        pt.time_reported = pt.time_symptom_onset + p.time_symptom_onset_to_detect; // Note double -> int conversion
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
    o.last_time_completed = 0;
    o.origin_time = 0;
}

void RunSimulationTime(const run_params &p, 
                       const int min_time, 
                       const int max_time, 
                       const int n_detections, 
                       std::vector<int> &t_detects, 
                       std::vector<int> &pop_size, // Population size at integer time t, relative to the index case
                       outbreak &o, 
                       gsl_rng *rgen)
{
    // This is the main simulation code
    int t = 0;                              // Time is an integer here ....
    MakeNewCase(p, -1, t_detects, o, rgen); // Make index case
    std::vector<int> t_detects_relative;    // Will store detection times relative to the first detection
    int index = -1;
    int zeros = 0; // Keep track of the number of generations with no new cases
    while (t < 1000)
    {
        index++;
        int added = 0;
        pop_size.push_back(0); // Records the number of people infected at exactly time t
        for (unsigned int i = 0; i < o.individuals.size(); i++)
        { 
            if (int(o.individuals[i].time_infected) == t) // Was Integer / double comparison
            {   // The population size stores the number of cases infected on this day.
                // pop_size will later be processed to get a total number of live cases, pop_sum.
                pop_size[index]++; 
            }
            if (int(o.individuals[i].time_symptom_onset) == t) // Was Integer / double comparison
            { // Generate infections for cases which become symptomatic on this day
                zeros = 0; // Reset zeros
                int infect = 0; // How many new cases generated by this individual
                if (added < p.add_limit)
                { // We have a limit here on how many cases to generate each day.
                    // The limit can be effectively removed with the --no_limit flag.
                    // Having a limit makes the code run faster for outbreaks with larger R0.  The default limit of 1000 ensures that some cases
                    // will almost certainly be detected each day, therefore altering the stats.
                    // Some thought needs to go into what happens if there are lots of detections each day in the data: Probably that isn't a good application of this code.
                    infect = gsl_ran_poisson(rgen, p.r0);
                    added += infect;
                }
                for (int j=0;j<infect;j++) {  //Add infect new cases of infection with parent index i
                    MakeNewCase (p,i,t_detects,o,rgen); /* These should all be infected at new time >=t */
                }
            }
        }
        // cout << "Number detected is now " << t_detects.size() << "\n";
        zeros++;
        o.last_time_completed = t;

        // Calculate relative detection times
        MakeRelativeTime(t_detects, t_detects_relative, o);

        // Check whether to terminate the simulation
        std::sort(t_detects_relative.begin(), t_detects_relative.end());
        int term = CheckTermination(p, t, zeros, min_time, max_time, n_detections, t_detects_relative, o);

        if (term == 1)
        {
            break;
        }
        /* What's consistent at this point (time t)
        we've counted all the individuals infected <= t

        */

        t++;
    }
    // Ensure detection times are sorted (Should be from l 139)
    std::sort(t_detects_relative.begin(), t_detects_relative.end());
    t_detects = t_detects_relative; // N.B. At the end of the simulation, the std::vector of detections is converted to relative time.
}

void MakeRelativeTime(const std::vector<int> &t_detects, 
                      std::vector<int> &t_detects_relative, 
                      outbreak &o)
{
    t_detects_relative = t_detects; // Copy the list of detection times
    int min = 1000;
    // Find min(min(t_detects), 1000)
    for (unsigned int i = 0; i < t_detects_relative.size(); i++)
    {
        if (t_detects_relative[i] < min)
        {
            min = t_detects_relative[i];
        }
    }
    // Subtract min fom all detection times
    for (unsigned int i = 0; i < t_detects_relative.size(); i++)
    {
        t_detects_relative[i] = t_detects_relative[i] - min;
    }
    // Fix the time of the index infection
    o.origin_time = -min;
}

int CheckTermination(const run_params &p, 
                     const int t, 
                     const int zeros, 
                     const int min_time, 
                     const int max_time, 
                     const int n_detections, 
                     const std::vector<int> &t_detects_relative, 
                     outbreak &o)
{
    // Terminate simulation if there have been more than infection_length+1 days with no individuals developing symptoms on that day
    // With Flu parameters this is probably OK, but it will require time to symptoms and time to infection to be almost always less than infection_length
    if (zeros > p.infection_length + 1)
    {
        if (p.verb == 1)
        {
            std::cout << "Terminate: Outbreak has died out\n";
        }
        // cout << "Zero out\n";
        return 1;
    }
    int crit_time = max_time; // Max time is the time of the last detection in the data.

    if (t_detects_relative.size() > n_detections)
    {
        crit_time = t_detects_relative[n_detections]; // First (relative) time at which the simulation has too many detections.
    }
    if (o.time_first_detect != -1 && t - o.time_first_detect >= crit_time) // ensure we have reached the time of the extra detection,
    {
        // What this is about: We have now calculated sufficiently many days worth of outbreak after the last detection that we can stop the simulation
        // We now have sufficient data to do anything we need to do subsequently with the data
        if (p.verb == 1)
        {
            std::cout << "Terminate: End of time after last detection\n";
        }
        return 1;
    }

    // Test whether the number of detections before the first measured time-point is greater than the number of detections in the dataset
    int det_early = 0;
    for (unsigned int i = 0; i < t_detects_relative.size(); i++)
    {
        if (t_detects_relative[i] <= min_time)
        { // Number of detections happening before the first evaluation point
            det_early++;
        }
    }

    if (det_early > n_detections)
    { // This number exceeds the total number of detections: Easy flag for non-compatibility with the data
        if (p.verb == 1)
        {
            std::cout << "Terminate: Too many detections too early\n";
        }
        return 1;
    }
    else
    { // NB There may be ways to increase efficiency by terminating things earlier but the above criteria are sufficient
        return 0;
    }
    
    
}

// This is weird - if t_detects.size()>0 and t_detects[0] !=0 then it inserts (0,0) into the list of detections
// However, think that for all inputs here t_detects[0] == 0.
void ConstructSummaryData(const std::vector<int> &t_detects, std::vector<detect> &sim_data)
{
    // Produces counts describing detections in the simulation in a format that matches the input data
    // Assume t_detects is already sorted
    int index = 0;
    int count = 0;
    //std::sort(t_detects.begin(), t_detects.end());
    for (unsigned int i = 0; i < t_detects.size(); i++)
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

void EvaluateOutbreak(const run_params &p, int exclude, const int r0val, std::vector<int> &t_detects, std::vector<int> &timepoints, std::vector<int> &pop_sum, std::vector<detect> &detections, outbreak &o, std::vector<std::vector<output>> &results)
{
    for (unsigned int i = 0; i < timepoints.size(); i++)
    {                               // Cycle through time points at which we are evaluating the simulation
        results[i][r0val].tested++; // Count of number of simulations
    }
    if (exclude == 0) // currently always the case
    {
        // Check the simulation against the detection data
        std::vector<detect> sim_data;
        ConstructSummaryData(t_detects, sim_data);
        if (p.verb == 1)
        {
            OutputSummaryData(sim_data);
        }
        /* Test for consistency at each observed timepoint */
        for (unsigned int i = 0; i < timepoints.size(); i++)
        {   
            int accept = 1;
            for (unsigned int j = 0; j < sim_data.size(); j++)
            {
                /* This looks suspect - could be a detection in the simulation after timepoint[i] */
                /*
                if (j >= detections.size())
                { // May be more detections in the simulation than in the data
                    accept = 0;
                    break;
                }
                else if (sim_data[j].day <= timepoints[i])
                {
                    if (sim_data[j].day != detections[j].day || sim_data[j].cases != detections[j].cases)
                    { // Mismatch with data
                        accept = 0;
                        break;
                    }
                }
                */
               // Restructure as -
                if (sim_data[j].day <= timepoints[i])
                {
                    if (j>=detections.size() || sim_data[j].day != detections[j].day || sim_data[j].cases != detections[j].cases)
                    { // Mismatch with data
                        accept = 0;
                        break;
                    }
                }              

            }
            if (accept == 1)
            { // Simulation fits the data
                if (p.verb == 1)
                {
                    std::cout << "Accepted at time point " << timepoints[i] << "\n";
                    std::cout << "Size of pop_sum " << pop_sum.size() << "\n";
                    int time = o.time_first_detect + timepoints[i];
                    std::cout << "Time " << time << "\n";
                    if (pop_sum.size() > time)
                    { // Current number of infections
                        std::cout << "Size " << pop_sum[time] << "\n";
                    }
                }
                results[i][r0val].accepted++;                           // Record accptance
                results[i][r0val].origin_time.push_back(o.origin_time); // Origin time of outbreak relative to detection
                // Find current size
                int time = o.time_first_detect + timepoints[i];
                if (pop_sum.size() > time)
                { // Store current number of "active" infections at this time
                    results[i][r0val].current_size.push_back(pop_sum[time]);
                }
                else
                {
                    results[i][r0val].current_size.push_back(0);
                    results[i][r0val].dead++;
                }
            }
        }
    }
    if (p.verb == 1)
    {
        std::cout << "\n";
    }
}

void MakePopulationSize(const run_params &p, const std::vector<int> &pop_size, std::vector<int> &pop_sum)
{
/*  Take input pop_size array. For every element index i, add pop_size[i] to pop_sum[i..i+p.infection_length-1] inclusive.
    Used to calculate the number of infections at each time, assuming that infections last from day x to day x+p.infection_length-1
*/
    for (unsigned int i = 0; i < pop_sum.size(); i++)
    {
        pop_sum[i] = 0;
    }
    for (unsigned int i = 0; i < pop_sum.size(); i++)
    {
        if (pop_size[i] > 0)
        {
            for (int j = 0; j < p.infection_length; j++)
            {
                if (i + j < pop_sum.size())
                {
                    pop_sum[i + j] = pop_sum[i + j] + pop_size[i];
                }
            }
        }
    }
}

void CalculateAcceptance(const run_params &p, const int i, const std::vector<std::vector<output>> &results, std::vector<double> &acceptance)
{
    double tot = 0;
    for (double r0 = 0.1; r0 <= 4.01; r0 = r0 + 0.1)
    {
        int r0val = floor((r0 + 0.001) * 10);
        double acc = (results[i][r0val].accepted + 0.) / (results[i][r0val].tested + 0.);
        tot = tot + acc;
        acceptance.push_back(acc);
    }
    for (unsigned int j = 0; j < acceptance.size(); j++)
    {
        acceptance[j] = acceptance[j] / tot;
    }
}
