#include "basicmodel.h"
#include "io.h"
#include "pseudorandom.h"
#include "utilities.h"
#include <algorithm>
#include <string>

void ConstructResults(vector<int>& timepoints, vector< vector<output> >& results) {
    for (unsigned int i=0;i<timepoints.size();i++) {
        vector<output> out;
        for (int i=0;i<51;i++) {
            output oo;
            oo.tested=0;
            oo.accepted=0;
            oo.dead=0;
            out.push_back(oo);
        }
        results.push_back(out);
    }
}

void MakeIndexCase (run_params& p, int& t, vector<int>& t_detects, outbreak& o, gsl_rng *rgen) {
    pat pt;
    pt.time_i=0;
    //pt.infected_by=-1;
    pt.time_s=pt.time_i+floor(gsl_ran_weibull(rgen,p.incubation_b,p.incubation_a)+0.5);
    o.last_time_completed=pt.time_s;
    pt.detected=gsl_ran_bernoulli(rgen,p.first_detect);
    if (pt.detected==1) {  //Could make more complicated by adding in a distribution for p.time_r?
        pt.time_r=pt.time_s+p.symptom_to_detect;
        o.first_detect=pt.time_r;
        t_detects.push_back(pt.time_r);
    } else {
        pt.time_r=-1;
    }
    o.indiv.push_back(pt);
}

void MakeIndexCaseFaster (run_params& p, int& t, const int& N, vector<int>& t_detects, const int& pp, long long& r, const long long& r_orig, vector<int>& new_incubate, vector<int>& new_detect, outbreak& o, gsl_rng *rgen) {
    pat pt;
    pt.time_i=0;
    //pt.infected_by=-1;
    int digit=GetRandomDigit(r_orig,pp,N,r);
    pt.time_s=pt.time_i+new_incubate[digit];
    o.last_time_completed=pt.time_s;
    pt.detected=gsl_ran_bernoulli(rgen,p.first_detect);
    if (pt.detected==1) {  //Could make more complicated by adding in a distribution for this parameter?
        pt.time_r=pt.time_s+p.symptom_to_detect;
        o.first_detect=pt.time_r;
        t_detects.push_back(pt.time_r);
    } else {
        pt.time_r=-1;
    }
    o.indiv.push_back(pt);
}


void MakeNewCase (run_params& p, int by, vector<int>& t_detects, outbreak& o, gsl_rng *rgen) {
    //Uses random numbers
    pat pt;
    //Time of infection
    pt.time_i=o.indiv[by].time_s+floor(gsl_ran_weibull(rgen,p.infection_b,p.infection_a)+0.5);
    pt.time_s=pt.time_i+floor(gsl_ran_weibull(rgen,p.incubation_b,p.incubation_a)+0.5);
    pt.detected=gsl_ran_bernoulli(rgen,p.detect);
    //cout << "Times " << "From " << by << " "  << o.indiv[by].time_s << " " << pt.time_i << " " << pt.time_s << " Detected ";
    if (pt.detected==1) {
        pt.time_r=pt.time_s+p.symptom_to_detect;
        if (o.first_detect==-1||pt.time_r<o.first_detect) {
            o.first_detect=pt.time_r;
        }
        t_detects.push_back(pt.time_r);
    } else {
        pt.time_r=-1;
    }
    o.indiv.push_back(pt);
}

void MakeNewCaseFaster  (run_params& p, int by, vector<int>& t_detects, const int& N, const int& pp, long long& r, const long long& r_orig, vector<int>& new_infect, vector<int>& new_incubate, vector<int>& new_detect, outbreak& o) {
    //Uses pseudorandom numbers
    pat pt;
    //Time of infection
    int digit=GetRandomDigit(r_orig,pp,N,r);
    pt.time_i=o.indiv[by].time_s+new_infect[digit];
    //Time of symptom onset
    digit=GetRandomDigit(r_orig,pp,N,r);
    pt.time_s=pt.time_i+new_incubate[digit];
    digit=GetRandomDigit(r_orig,pp,N,r);
    pt.detected=new_detect[digit];
    if (pt.detected==1) {
        pt.time_r=pt.time_s+p.symptom_to_detect;
        if (o.first_detect==-1||pt.time_r<o.first_detect) {
            o.first_detect=pt.time_r;
        }
        t_detects.push_back(pt.time_r);
    } else {
        pt.time_r=-1;
    }
    o.indiv.push_back(pt);
}

void SetupOutbreak(outbreak& o) {
    o.first_detect=-1;
    o.total_detections=0;
    o.last_time_completed=0;
    o.origin_time=0;
}

void RunSimulationTime (run_params& p, int min_time, int max_time, int& exclude, const int& n_detections, const int& N, const int& pp, long long& r, const long long& r_orig, vector<int>& new_infect, vector<int>& new_incubate, vector<int>& new_detect, vector<int>& new_number, vector<int>& t_detects, vector<int>& pop_size, outbreak& o, gsl_rng *rgen) {
    //This is the main simulation code
    int t=0;
    if (p.run_fast==1) {
        MakeIndexCaseFaster (p,t,N,t_detects,pp,r,r_orig,new_incubate,new_detect,o,rgen);
    } else {
        MakeIndexCase (p,t,t_detects,o,rgen);
    }
    vector<int> t_detects_relative;  //Will store detection times relative to the first detection being on day 0
    int index=-1;
    int zeros=0; //Keep track of the number of generations with no new cases
    while (t<1000) {
        index++;
        int added=0;
        int infect=0;
        pop_size.push_back(0);  //Records the number of people infected at exactly time t
        //cout << "T " << t << " " << pop_size.size() << " Zeros " << zeros << " Size " << o.indiv.size() << " Detected " << t_detects.size() << " First " << o.first_detect << "\n";
        for (unsigned int i=0;i<o.indiv.size();i++) {  //Keep going as there could be cases with instant symptom onset
            if (o.indiv[i].time_i==t) { //The population size stores the number of cases infected on this day.
                //pop_size will later be processed to get a total number of live cases, pop_sum.
                pop_size[index]++;
            }
            if (o.indiv[i].time_s==t) { //Generate infections for cases which become symptomatic on this day
                zeros=0;
                if (added<p.add_limit) { //We have a limit here on how many cases to generate each day.
                    //The limit can be effectively removed with the --no_limit flag.
                    //Having a limit makes the code run faster for outbreaks with larger R0.  The default limit of 1000 ensures that some cases
                    //will almost certainly be detected each day, therefore altering the stats.
                    //Some thought needs to go into what happens if there are lots of detections each day in the data: Probably that isn't a good application of this code.
                    if (p.run_fast==1) {
                        int digit=GetRandomDigit(r_orig,pp,N,r);
                        infect=new_number[digit];
                    } else {
                        infect=gsl_ran_poisson(rgen,p.r0);
                    }
                    added=added+infect;
                } else {
                    infect=0; //At the limit, just don't generate any more cases.  I don't think there is anything special about the first cases with p.time_s=t.
                }
                //cout << "Individual " << i << " Infects " << infect << " individuals\n";
                for (int j=0;j<infect;j++) {  //Add new cases of infection
                    //o.indiv[i].infects.push_back(o.indiv.size());  //Conceivably don't need this line.
                    if (p.run_fast==1) {
                        MakeNewCaseFaster(p,i,t_detects,N,pp,r,r_orig,new_infect,new_incubate,new_detect,o);
                    } else {
                        MakeNewCase (p,i,t_detects,o,rgen);
                    }
                }
            }
        }
        //cout << "Number detected is now " << t_detects.size() << "\n";
        zeros++;
        o.last_time_completed=t;
                
        //Calculate relative detection times
        MakeRelativeTime (t_detects,t_detects_relative,o);
        
        //Check whether to terminate the simulation
        sort(t_detects_relative.begin(),t_detects_relative.end());
        int term=CheckTermination(p,t,zeros,min_time,max_time,n_detections,t_detects_relative,o);
        if (term==1) {
            break;
        }
        t++;
    }
    t_detects=t_detects_relative; //N.B. At the end of the simulation, the vector of detections is converted to relative time.
}

void MakeRelativeTime (vector<int>& t_detects, vector<int>& t_detects_relative, outbreak& o) {
    t_detects_relative=t_detects;
    int min=1000;
    for (unsigned int i=0;i<t_detects_relative.size();i++) {
        if (t_detects_relative[i]<min) {
            min=t_detects_relative[i];
        }
    }
    for (unsigned int i=0;i<t_detects_relative.size();i++) {
        t_detects_relative[i]=t_detects_relative[i]-min;
    }
    o.origin_time=-min;
}

int CheckTermination (run_params& p, int& t, int& zeros, int& min_time, int& max_time, const int& n_detections, vector<int>& t_detects_relative, outbreak& o) {
    if (zeros>p.infection_length+1) {
        if (p.verb==1) {
            cout << "Terminate: Outbreak has died out\n";
        }
        //cout << "Zero out\n";
        return 1;
    } else {
        int crit_time=max_time; //Max time is the time of the last detection in the data.
        if (t_detects_relative.size()>n_detections) {
            crit_time=t_detects_relative[n_detections];
        }
        if (o.first_detect!=-1&&t-o.first_detect>=crit_time) {
            //What this is about: We have now calculated sufficiently many days worth of outbreak after the last detection that we can stop the simulation
            //We now have sufficient data to do anything we need to do subsequently with the data
            if (p.verb==1) {
                cout << "Terminate: End of time after last detection\n";
            }
            return 1;
        } else { //Too many detections prior to the first evaluation point?
            int det_early=0;
            for (unsigned int i=0;i<t_detects_relative.size();i++) {
                if (t_detects_relative[i]<=min_time) { //Number of detections happening before the first evaluation point
                    det_early++;
                }
            }
            if (det_early>n_detections) { //This number exceeds the total number of detections: Easy flag for non-compatibility with the data
                if (p.verb==1) {
                    cout << "Terminate: Too many detections too early\n";
                }
                return 1;
            } else {  //NB There may be ways to increase efficiency by terminating things earlier but the above criteria are sufficient
                return 0;
            }
        }
    }
}

void ConstructSummaryData (vector<int>& t_detects, vector<detect>& sim_data) {
    //Produces counts describing detections in the simulation in a format that matches the input data
    int index=0;
    int count=0;
    sort(t_detects.begin(),t_detects.end());
    for (unsigned int i=0;i<t_detects.size();i++) {
        if (t_detects[i]==index) {
            count++;
        } else {
            detect d;
            d.day=index;
            d.cases=count;
            sim_data.push_back(d);
            index=t_detects[i];
            count=1;
        }
    }
    detect d;
    d.day=index;
    d.cases=count;
    sim_data.push_back(d);
}

void EvaluateOutbreak (run_params& p, int& exclude, int& r0val, vector<int>& t_detects, vector<int>& timepoints, vector<int>& pop_sum, vector<detect>& detections, outbreak& o, vector< vector<output> >& results) {
    for (unsigned int i=0;i<timepoints.size();i++) { //Cycle through time points at which we are evaluating the simulation
        results[i][r0val].tested++; //Count of number of simulations
    }
    if (exclude==0) {
        //Check the simulation against the detection data
        vector<detect> sim_data;
        ConstructSummaryData (t_detects, sim_data);
        if (p.verb==1) {
            OutputSummaryData (sim_data);
        }
        for (unsigned int i=0;i<timepoints.size();i++) { //Cycle through time points at which we are evaluating the simulation
            //These represent different periods of time after which no more cases have been observed
            int accept=1;
            for (unsigned int j=0;j<sim_data.size();j++) {
                if (sim_data[j].day<=timepoints[i]) {
                    if (j>=detections.size()) { //May be more detections in the simulation than in the data
                        accept=0;
                        break;
                    } else if (sim_data[j].day!=detections[j].day||sim_data[j].cases!=detections[j].cases) { //Mismatch with data
                        accept=0;
                        break;
                    }
                }
            }
            if (accept==1) { //Simulation fits the data
                if (p.verb==1) {
                    cout << "Accepted at time point " << timepoints[i] << "\n";
                    cout << "Size of pop_sum " << pop_sum.size() << "\n";
                    int time=o.first_detect+timepoints[i];
                    cout << "Time " << time << "\n";
                    if (pop_sum.size()>time) { //Current number of infections
                        cout << "Size " << pop_sum[time] << "\n";
                    }
                }
                results[i][r0val].accepted++;  //Record accptance
                results[i][r0val].origin_time.push_back(o.origin_time); //Origin time of outbreak relative to detection
                //Find current size
                int time=o.first_detect+timepoints[i];
                if (pop_sum.size()>time) { //Current number of infections
                    results[i][r0val].current_size.push_back(pop_sum[time]);
                } else {
                    results[i][r0val].current_size.push_back(0);
                    results[i][r0val].dead++;
                }
            }
        }
    }
    if (p.verb==1) {
        cout << "\n";
    }
}


void MakePopulationSize (run_params& p, vector<int>& pop_size, vector<int>& pop_sum) {
    for (unsigned int i=0;i<pop_sum.size();i++) {
        pop_sum[i]=0;
    }
    for (unsigned int i=0;i<pop_sum.size();i++) {
        if (pop_size[i]>0) {
            for (int j=0;j<p.infection_length;j++) {
                if (i+j<pop_sum.size()) {
                    pop_sum[i+j]=pop_sum[i+j]+pop_size[i];
                }
            }
        }
    }
}

void CalculateAcceptance (run_params& p, int i, const vector< vector<output> >& results, vector<double>& acceptance) {
    double tot=0;
    for (p.r0=0.1;p.r0<=4.01;p.r0=p.r0+0.1) {
        int r0val=floor((p.r0+0.001)*10);
        double acc=(results[i][r0val].accepted+0.)/(results[i][r0val].tested+0.);
        tot=tot+acc;
        acceptance.push_back(acc);
    }
    for (unsigned int j=0;j<acceptance.size();j++) {
        acceptance[j]=acceptance[j]/tot;
    }
}
