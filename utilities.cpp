#include "basicmodel.h"
#include "io.h"
#include "utilities.h"
#include <algorithm>
#include <string>

void ConstructResults(const run_params &p, const std::vector<int> &timepoints, std::vector<std::vector<output> > &results)
{

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

inline void MakeNewCase(const run_params &p, const int by, std::vector<int> &t_detects, std::vector<int> &t_detects_after, outbreak &o, gsl_rng *rgen)
{


    patient pt;
    // Time of infection
    if (by == -1) // Flag for index case
    {
        pt.time_infected = 0;
    }
    else
    {
      pt.time_infected = o.individuals[by].time_symptom_onset + static_cast<int>(floor(gsl_ran_weibull(rgen, p.infection_b, p.infection_a)+0.5)); 
    }
    pt.time_symptom_onset = pt.time_infected + static_cast<int>(floor(gsl_ran_weibull(rgen, p.incubation_b, p.incubation_a)+0.5));
    pt.infection_length = gsl_rng_uniform_int(rgen, p.infection_length_max - p.infection_length_min + 1) + p.infection_length_min;

    if (by == -1) // Index case
    {
        double r= gsl_rng_uniform(rgen);
        pt.detected = (r<=p.probability_first_detect);
        pt.detected_after = (r<= p.probability_detect_enhanced); // For some parameters, it is possible (but unlikely)
	                                                         // for another case to be detected before this one could be detected.
    }
    else
    {
        // replace sample
        double r = gsl_rng_uniform(rgen);
        pt.detected = (r<=p.probability_detect);
        pt.detected_after = (r<= p.probability_detect_enhanced); // Below code assumes probability_detected <= probability_detected_after
    }
    if ((pt.detected == 1) || (pt.detected_after==1))
    {
        if (pt.detected==1) {
	  int time_symptom_onset_to_detect = gsl_rng_uniform_int(rgen, p.time_symptom_onset_to_detect_max - p.time_symptom_onset_to_detect_min + 1) + p.time_symptom_onset_to_detect_min;
	  pt.time_reported = pt.time_symptom_onset + time_symptom_onset_to_detect; 
            if (o.time_first_detect == -1 || pt.time_reported < o.time_first_detect)
            {
                o.time_first_detect = pt.time_reported;
		/*
                for(int t : t_detects_after) {
                    if(t>=o.time_first_detect) {
                        t_detects.push_back(t);
                    }
                } 
		*/
                auto should_remove = [o,p](int t) { return t>=o.time_first_detect+p.time_before_enhanced_detection; }; 

                auto new_end = std::remove_if(t_detects_after.begin(), t_detects_after.end(), should_remove);

		t_detects.insert(t_detects.end(), new_end, t_detects_after.end());
                t_detects_after.erase(new_end, t_detects_after.end());
            }
            t_detects.push_back(pt.time_reported);
        } else {
	  int time_symptom_onset_to_detect = gsl_rng_uniform_int(rgen, p.time_symptom_onset_to_detect_max - p.time_symptom_onset_to_detect_min + 1) + p.time_symptom_onset_to_detect_min;

	    pt.time_reported = pt.time_symptom_onset + time_symptom_onset_to_detect; 
            if (o.time_first_detect != -1 && pt.time_reported >= o.time_first_detect + p.time_before_enhanced_detection)
            {
                t_detects.push_back(pt.time_reported);

            } else {
                t_detects_after.push_back(pt.time_reported);
            }
        }
    }
    else
    {
        pt.time_reported = -1;
    }
    o.individuals.push_back(pt);
}

void SetupOutbreak(outbreak &o)
{
    o.time_first_detect = -1; 
    o.last_time_simulated = 0;
}

void RunSimulationTime(const run_params &p, 
                       const double r0,
                       const int first_evaluation_time, 
                       const int last_evaluation_time, 
                       const unsigned long n_data_detections, 
                       const int extreme_infection_time,
                       std::vector<int> &t_detects_relative,
                       std::vector<int> &number_new_symptomatic, 
                       outbreak &o, 
                       gsl_rng *rgen)
{


    int t = 0;
    std::vector<int> t_detects;
    std::vector<int> t_detects_after;
    MakeNewCase(p, -1, t_detects, t_detects_after, o, rgen); // Make index case
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
		  MakeNewCase (p,i,t_detects,t_detects_after, o,rgen); /* These should all be infected at new time >=t */
                }
            }
        }
        zeros++;
        o.last_time_simulated = t;

        // Check whether to terminate the simulation (early)
        std::sort(t_detects.begin(), t_detects.end());
        int term = CheckTermination(p, t, zeros, extreme_infection_time, first_evaluation_time, last_evaluation_time, n_data_detections, t_detects, o);

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
}

int CheckTermination(const run_params &p, 
                     const int t, 
                     const int zeros, 
                     const int extreme_infection_time,
                     const int first_evaluation_time,
                     const int last_evaluation_time,
                     const unsigned long n_detections, 
                     const std::vector<int> &t_detects, 
                     const outbreak &o)
{
    /** Each simulation proceeded up to the either 
     * a) the outbreak had certainly died out or
     * b) it was clear that the simulation would either be accepted or rejected at all observation times.
     *    Note that rejection at one observation time implies rejection at all later observation times.
     *    We also need to have completely simulated up to the last observation time that might be accepted (for other statistics).
     * 
     * Here, for b) we terminated simulations when either:
     * i) The time t has reached or gone beyond the end of the day t_detect[0]+last_observation_time. (Everything is known.)
     *    Then the first detected case must be before t, and there is no possibility of any additional detections before t.
     * ii) There are more detected cases than in the dataset, and that t has reached the end of the day in which the
     *     first extra case was detected. (Know any observation timepoint after t will be rejected, and any before is completely simulated)
     * iii) There too many detected cases before the earliest time that the first observation timepoint could occur. (All rejected.)   
     * 
     * More optimal strategies are likely possible, especially in cases where the experimental data has more than one detected case.
     */

    // a) Terminate simulation if there have been more than extreme_infection_time days with no individuals developing symptoms on that day (zeros is the number of days +1 since an individual last developed symptoms)
    if (zeros > extreme_infection_time ) 
    {
        if (p.verb == 1)
        {
            std::cout << "Terminate: Outbreak has died out at time " << t << "zeros: " << zeros <<"\n";
        }
        // cout << "Zero out\n";
        return 1;
    }

    if ((o.time_first_detect != -1) && (t - o.time_first_detect >= last_evaluation_time)) 
    {
        // b) i) This stops the simulation once the time has reached the last observation. 
        // (At this point t >= time_first_detect, so first detection is known)
        // Matching will not consider any cases after this time
        if (p.verb == 1)
        {
            std::cout << "Terminate: End of time " << t << " after last detection\n";
        }
        return 1;
    }  



    if (t_detects.size() > n_detections)
    {
        if ((o.time_first_detect != -1) && (t  >= t_detects[n_detections])) 
        {
            // b) ii) At this point we know that the first detection is at t_detects[0],
            // all timepoints on or after t_detects[n_detections]
            // will fail, and all those before are completely simulated.
            if(p.verb==1)
            {
                std::cout << "Terminate: Too many detections "<< t_detects.size() << " vs " << n_detections << " before time " << t <<"\n";
            }
            return 1;
        }
    }

    // b) iii) There too many cases before the earliest time that the first observation timepoint could occur.
    // In this case, we know that all the matching comparisons with the data will fail.
    // Simulation could generate a new detected case at t_new_detect in (t, t_detect[0]]
    // which might be match the data for [t_new_detect, t_new_detect+first_evaluation_time]
    // Know that the first detection occurs after t'=min(t, t_detect[0])
    // so count detections in range [min(t, t_detect[0]), min(t, t_detect[0])+first_evaluation_time]
    // If there are too many, then all timepoints will be rejected 

    // Could potentially improve on this knowing time_symptom_onset_to_detect_min, 
    // Simulation could generate a new detected case in (t+p.time_symptom_onset_to_detect_min, t_detect[0]]
    // so set t' = min (t+p.time_symptom_onset_to_detect_min+1, t_detect[0]) - note also sharper bound on t, as we have finished simulating all cases which became symptomatic at time t
    // count detections in range [t', t'+first_evaluation_time] -> very little improvement in simulation time
    // This bound on the number of cases could also be tightened using the data - here it is just the total number of detections.


    if((o.time_first_detect != -1)) {
        // Test whether the number of detections on or before the first evaluation time-point is greater than the number of detections in the dataset - if so, this will be the case for all later ones.
        unsigned long det_early = 0;

        // Now assume sorted - break if time greater than detection time
        for (unsigned long i = 0; (i < t_detects.size()) && (t_detects[i]<=(std::min(o.time_first_detect/*+1+p.time_symptom_onset_to_detect*/, t) + first_evaluation_time)); i++)
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



    if(t_detects.empty()) {
        return;
    }
    // Produces counts describing detections in the simulation in a format that matches the input data
    // Assume t_detects is already sorted
    int index = 0;
    int count = 0;
    for (unsigned long i = 0; i < t_detects.size(); i++)
    {
        if (t_detects[i] == index)
        {
            count++;
        }
        else
        {
            if(count>0)
            {
                detect d;
                d.day = index;
                d.cases = count;
                sim_data.push_back(d);
            }
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

void EvaluateOutbreak(const run_params &p, const unsigned long r0val, std::vector<int> &t_detects_relative, std::vector<int> &timepoints, std::vector<int> &total_active_infected, std::vector<detect> &detections, outbreak &o, std::vector<std::vector<output> > &results)
{


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
        if(o.time_first_detect<0) {
            std::cout << "Negative first detection time " << o.time_first_detect << " " << n_timepoints_accepted << " " << sim_data.size() << "\n";
        }
        // Find the total number of currently infected individuals at the observation timepoint
        int time = o.time_first_detect + timepoints[i];
        if (total_active_infected.size() > time)
        { // Store current number of "active" infections at this time
            results[i][r0val].current_size.push_back(total_active_infected[time]);
        }
        else
        {
            results[i][r0val].current_size.push_back(0);
            results[i][r0val].dead++; // No currently infected individuals at the observation timepoint
        }
}
    if (p.verb == 1)
    {
        std::cout << "\n";
    }
}

void MakePopulationSize(const run_params &p, const outbreak &o, std::vector<int> &total_active_infected)
{


    for (unsigned long i = 0; i < total_active_infected.size(); i++)
    {
        total_active_infected[i] = 0;
    }
    for(unsigned long i = 0; i < o.individuals.size(); i++)
      {
	for (unsigned long j=o.individuals[i].time_symptom_onset; j < o.individuals[i].time_symptom_onset + o.individuals[i].infection_length && j < total_active_infected.size(); j++)
	  {
	    total_active_infected[j]++;
	  }
      }
}
