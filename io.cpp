#include "basicmodel.h"
#include "io.h"
#include "utilities.h"
#include <string>
#include <algorithm>


void GetOptions (run_params& p, int argc, const char **argv) {
  std::string p_switch;
    p.input_prefix="Data/";
    p.output_prefix="";
    p.incubation_a=7.402580682098853; // Too many SF
    p.incubation_b=1.7375081458523993;
    p.infection_a=1.03138989929643;
    p.infection_b=1.0025194828961315;
    p.max_infections=500000;
    p.add_limit=1000;
    p.probability_detect=0.1;
    p.probability_first_detect=-1; //Probability of first detection.  By default set to p.detect later.
    p.time_symptom_onset_to_detect=18;
    p.replicas=1000000;
    p.infection_length=7;
    p.max_R0=4.0;
    p.more_stats=0;
    p.verb=0;
    p.seed=1234;
    p.max_simulation_time=1000;
    int x=1;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
        if (p_switch.compare("--max_infections")==0) {
            x++;
            p.max_infections=atoi(argv[x]);
        } else if (p_switch.compare("--replicas")==0) {
            x++;
            p.replicas=atoi(argv[x]);
        } else if (p_switch.compare("--detect")==0) {
            x++;
            p.probability_detect=atof(argv[x]);
        } else if (p_switch.compare("--first_detect")==0) {
            x++;
            p.probability_first_detect=atof(argv[x]);
        } else if (p_switch.compare("--verb")==0) {
            x++;
            p.verb=atoi(argv[x]);
        } else if (p_switch.compare("--max_R0")==0) {
            x++;
            p.max_R0=atof(argv[x]);
        } else if (p_switch.compare("--more_stats")==0) {
            x++;
            p.more_stats=atoi(argv[x]);
        } else if (p_switch.compare("--no_limit")==0) {  //Remove the limit on the population size
            x++;
            p.add_limit=pow(10,9);
        } else if (p_switch.compare("--output_prefix")==0) {  
            x++;
            p.output_prefix=argv[x];
        } else if (p_switch.compare("--input_prefix")==0) { 
            x++;
            p.input_prefix=argv[x];
        } else {
		  std::cout << "Incorrect usage " << argv[x] << "\n ";
		  exit(1);
		}
		p_switch.clear();
		x++;
	}
    p.R0_vals.clear();
    unsigned long  top=floor((p.max_R0*10)+0.5);
    for (unsigned long r0val=1; r0val<=top ;r0val++)
    {
        p.R0_vals.push_back(r0val*0.1);
    }
	if (p.probability_first_detect==-1) 
    {
        p.probability_first_detect=p.probability_detect;
    }
}

void ImportDetections(const run_params& p, int& n_detections, std::vector<detect>& detections)
{
    std::ifstream detect_file;
    detect_file.open(p.input_prefix+"Detections.dat");
    int day;
    int cases;
    for (unsigned long i=0;i<1000000;i++) {
        if (!(detect_file >> day)) break;
        if (!(detect_file >> cases)) break;
        if((day < 0) || (cases<=0))
        {
            std::cout << "Invalid data item: " << i << "day: " << day << "cases: " << cases << "\n";
            return;
        }
        detect d;
        d.day=day;
        d.cases=cases; // Note that we need this to always be bigger than zero, otherwise issues with comparison
        detections.push_back(d);
    }
    for(unsigned long i=0;i<detections.size();i++) {
        n_detections=n_detections+detections[i].cases;
    }
}

void ImportTimePoints(const run_params& p, int& min_time, int& max_time, std::vector<int>& timepoints)
{
    std::ifstream time_file;
    time_file.open(p.input_prefix+"Time_points.dat"); // Need these to be sorted in increasing order - so we sort them here. Require all to be positive (I think)
    int t;
    for (unsigned long i=0;i<1000000;i++)
    {
      if (!(time_file >> t)) break;
        if(t<0) {
            std::cout << "Invalid timepoint value t:" << t << "\n";
        }
        timepoints.push_back(t);
    }
    if(timepoints.size())
    {
        std::sort(timepoints.begin(), timepoints.end());
        max_time = timepoints[timepoints.size()-1];
        min_time = timepoints[0];
    } else {
        min_time = 0;
        max_time = 1000;
    }
}

void OutputSummaryData (const std::vector<detect>& sim_data) {
    //Matches the detections that were read in
    std::cout << "Sim data\n";
    for (unsigned long i=0;i<sim_data.size();i++) {
        std::cout << sim_data[i].day << " " << sim_data[i].cases << "\n";
    }
}

void OutputPopulationDetails (const run_params& p, const std::vector<int>& number_new_symptomatic, const std::vector<int>& total_active_infected, const std::vector<int> t_detects_relative, outbreak& o) {
    std::cout << "Population size\n";
    std::cout << number_new_symptomatic.size() << "\n";
    for (unsigned long i=0;i<number_new_symptomatic.size();i++) {
        std::cout << number_new_symptomatic[i] << " ";
    }
    std::cout << "\n";

    std::cout << "Population sum\n";
    std::cout << total_active_infected.size() << "\n";
    for (unsigned long i=0;i<total_active_infected.size();i++) {
        std::cout << total_active_infected[i] << " ";
    }
    std::cout << "\n";

    std::cout << "Number of cases " << o.individuals.size() << "\n";
    std::cout << "Number of detections " << t_detects_relative.size() << "\n";
    std::cout << "First detection " << o.time_first_detect << " of case symptomatic at time " << o.time_first_detect-p.time_symptom_onset_to_detect << "\n";
    std::cout << "T detects relative to first detection \n";
    for (unsigned long i=0;i<t_detects_relative.size();i++) {
        std::cout << t_detects_relative[i] << " ";
    }
    std::cout << "\n";

}

void OutputRawData (const int& r0val, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results) {
    for (unsigned long i=0;i<timepoints.size();i++) {
        std::cout << "Timepoint " << timepoints[i] << "\n";
        std::cout << "Total accepted " << results[i][r0val].accepted << "\n";
        std::cout << "Origin times\n";
        for (unsigned long j=0;j<results[i][r0val].origin_time.size();j++) {
            std::cout << results[i][r0val].origin_time[j] << " ";
        }
        std::cout << "\n";
        std::cout << "Population sizes\n";
        for (unsigned long j=0;j<results[i][r0val].current_size.size();j++) {
            std::cout << results[i][r0val].current_size[j] << " ";
        }
        std::cout << "\n";
    }
}

void OutputAcceptanceRates(const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results) {
    for (unsigned long r0val=0; r0val<p.R0_vals.size();r0val++) {
        std::ofstream acc_file;
        std::ostringstream convert;
        convert << r0val+1;
	    std::string temp=convert.str();
	    std::string name =p.output_prefix +  "Acceptance_rate"+temp+".dat";
        acc_file.open(name.c_str());
        for (unsigned long i=0;i<timepoints.size();i++) {
            double acc=static_cast<double>(results[i][r0val].accepted)/results[i][r0val].tested;
            acc_file << timepoints[i] << " " << acc << "\n";
        }
    }
}

void OutputOriginTimes (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results) {

    for (unsigned long r0val=0; r0val<p.R0_vals.size();r0val++) {
      std::ofstream init_file;
      std::ostringstream convert;
        convert << r0val+1;
	std::string temp=convert.str();
	std::string name =p.output_prefix +  "Origin_times"+temp+".dat";
        init_file.open(name.c_str());
        for (unsigned long i=0;i<timepoints.size();i++) {
            init_file << timepoints[i] << " ";
            for (unsigned long j=0;j<results[i][r0val].origin_time.size();j++) {
                init_file << results[i][r0val].origin_time[j] << " ";
            }
            init_file << "\n";
        }
        init_file.close();
    }
}

void OutputPopulationSizes (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results) {

    for (unsigned long r0val=0; r0val<p.R0_vals.size();r0val++) {

        std::ofstream size_file;
        std::ostringstream convert;
                convert << r0val+1;
	std::string temp=convert.str();
	std::string name = p.output_prefix + "Current_size"+temp+".dat";
        size_file.open(name.c_str());
        for (unsigned long i=0;i<timepoints.size();i++) {
            size_file << timepoints[i] << " ";
            for (unsigned long j=0;j<results[i][r0val].current_size.size();j++) {
                if (results[i][r0val].current_size[j]>0) {
                    size_file << results[i][r0val].current_size[j] << " ";
                }
            }
            size_file << "\n";
        }
        size_file.close();
    }
}

void OutputProbabilityEnded (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results) {

    for (unsigned long r0val=0; r0val<p.R0_vals.size();r0val++) {

      std::ofstream dead_file;
      std::ostringstream convert;
        convert << r0val+1;
	std::string temp=convert.str();
	std::string name = p.output_prefix + "End_prob"+temp+".dat";
        dead_file.open(name.c_str());
        for (unsigned long i=0;i<timepoints.size();i++) {
            double acc=(results[i][r0val].dead+0.)/(results[i][r0val].accepted+0.);
            dead_file << timepoints[i] << " " << acc << "\n";
        }
        dead_file.close();
    }
}

void OutputOutbreakDeathStatistics (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results) {
  std::ofstream pd_file;
    pd_file.open(p.output_prefix + "P_Outbreak_End.dat");
    for (unsigned long i=0;i<timepoints.size();i++) { //For each of the time points
        //Make acceptance std::vector for this timepoint
        std::vector<double> acceptance;
        CalculateAcceptance (p,i,results,acceptance);
        std::vector<double> pdead;
        for (unsigned long r0val=0; r0val<p.R0_vals.size();r0val++) {
            double pd=(results[i][r0val].dead+0.)/(results[i][r0val].accepted+0.);
            pdead.push_back(pd);
        }
        double prob=0;
        for (unsigned long j=0;j<acceptance.size();j++) {
            prob=prob+acceptance[j]*pdead[j];
        }
        pd_file << timepoints[i] << " " << prob << "\n";
    }
    pd_file.close();
}

void OutputOutbreakTimeStatistics (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results) {
    for (unsigned long i=0;i<timepoints.size();i++) { //For each of the time points
        //std::cout << "Time point " << timepoints[i] << "\n";
        std::vector<double> acceptance;
        CalculateAcceptance (p,i,results,acceptance);
        //Make a count of the origin dates for each R0, weighted by the acceptance rate for that R0
        std::vector< std::vector<double> > all_origins;
        for (unsigned long r0val=0; r0val<p.R0_vals.size();r0val++) {

            std::vector<double> origins;
            //std::cout << "R0 " << r0val << " Size " << results[i][r0val].origin_time.size() << "\n";
            for (unsigned long k=0;k<results[i][r0val].origin_time.size();k++) {
                if (results[i][r0val].origin_time[k]<=0) {
                    while (-results[i][r0val].origin_time[k]>=static_cast<int>(origins.size())) {
                        origins.push_back(0);
                    }
                    origins[-results[i][r0val].origin_time[k]]++;
                }
            }
            for (unsigned long j=0;j<origins.size();j++) {
                origins[j]=origins[j]*acceptance[r0val];
            }
            all_origins.push_back(origins);
        }
        //Compile all of these
        std::vector<double> combined_origins;
        for (unsigned long j=0;j<all_origins.size();j++) {
            for (unsigned long k=0;k<all_origins[j].size();k++) {
                while (k>=combined_origins.size()) {
                    combined_origins.push_back(0);
                }
                combined_origins[k]=combined_origins[k]+all_origins[j][k];
            }
        }
        
        double tot=0;
        for (unsigned long j=0;j<combined_origins.size();j++) {
            tot=tot+combined_origins[j];
        }
        for (unsigned long j=0;j<combined_origins.size();j++) {
            combined_origins[j]=combined_origins[j]/tot;
        }
	std::ofstream orig_file;
	std::ostringstream convert;
        convert << timepoints[i];
	std::string temp=convert.str();
	std::string name=p.output_prefix + "Origin_stats_"+temp+".dat";
        orig_file.open(name.c_str());
        for (unsigned long j=0;j<combined_origins.size();j++) {
            orig_file << j << " " << combined_origins[j] << "\n";
        }
        //std::cout << "\n";
        orig_file.close();
    }
}

void OutputOutbreakPopulationStatistics (const run_params& p, const std::vector<int>& timepoints, const std::vector< std::vector<output> >& results) {
    //Size of population conditional on not being zero
    for (unsigned long i=0;i<timepoints.size();i++) { //For each of the time points
        std::vector<double> acceptance;
        CalculateAcceptance (p,i,results,acceptance);
        //Make a count of the population sizes for each R0, weighted by acceptance rate
        std::vector< std::vector<double> > all_sizes;
        for (unsigned long r0val=0 ;r0val<p.R0_vals.size();r0val++) {
            std::vector<double> sizes;
            //std::cout << "R0 " << r0val << " Size " << results[i][r0val].current_size.size() << "\n";
            for (unsigned long k=0;k<results[i][r0val].current_size.size();k++) {
                while (results[i][r0val].current_size[k]>=sizes.size()) {
                    sizes.push_back(0);
                }
                if (results[i][r0val].current_size[k]>0) {
                    sizes[results[i][r0val].current_size[k]]++;
                }
            }
            for (unsigned long j=0;j<sizes.size();j++) {
                sizes[j]=sizes[j]*acceptance[r0val];
            }
            all_sizes.push_back(sizes);
        }
        //Compile all of these
        std::vector<double> combined_sizes;
        for (unsigned long j=0;j<all_sizes.size();j++) {
            for (unsigned long k=0;k<all_sizes[j].size();k++) {
                while (k>=combined_sizes.size()) {
                    combined_sizes.push_back(0);
                }
                combined_sizes[k]=combined_sizes[k]+all_sizes[j][k];
            }
        }
        double tot=0;
        for (unsigned long j=0;j<combined_sizes.size();j++) {
            tot=tot+combined_sizes[j];
        }
        for (unsigned long j=0;j<combined_sizes.size();j++) {
            combined_sizes[j]=combined_sizes[j]/tot;
        }
	std::ofstream size_file;
	std::ostringstream convert;
        convert << timepoints[i];
	std::string temp=convert.str();
	std::string name=p.output_prefix + "Sizes_time_"+temp+".dat";
        size_file.open(name.c_str());
        for (unsigned long j=0;j<combined_sizes.size();j++) {
            size_file << j << " " << combined_sizes[j] << "\n";
        }
        size_file.close();
    }
}
