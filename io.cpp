#include "basicmodel.h"
#include "io.h"
#include "utilities.h"
#include <string>

void GetOptions (run_params& p, int argc, const char **argv) {
  std::string p_switch;
    p.r0=1;
    p.incubation_a=7.402580682098853;
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
    p.more_stats=0;
    p.verb=0;
    p.seed=1234;
    int x=1;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--r0")==0) {
			x++;
			p.r0=atof(argv[x]);
        } else if (p_switch.compare("--max_infections")==0) {
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
        } else if (p_switch.compare("--more_stats")==0) {
            x++;
            p.more_stats=atoi(argv[x]);
        } else if (p_switch.compare("--no_limit")==0) {  //Remove the limit on the population size
            x++;
            p.add_limit=pow(10,9);

        } else {
		  std::cout << "Incorrect usage " << argv[x] << "\n ";
		  exit(1);
		}
		p_switch.clear();
		x++;
	}
	if (p.probability_first_detect==-1) {
        p.probability_first_detect=p.probability_detect;
    }
}

void ImportDetections(int& n_detections, std::vector<detect>& detections) {
    std::ifstream detect_file;
    detect_file.open("Data/Detections.dat");
    int day;
    int cases;
    for (int i=0;i<1000000;i++) {
        if (!(detect_file >> day)) break;
        if (!(detect_file >> cases)) break;
        detect d;
        d.day=day;
        d.cases=cases;
        detections.push_back(d);
    }
    for(int i=0;i<detections.size();i++) {
        n_detections=n_detections+detections[i].cases;
    }
}

void ImportTimePoints(int& min_time, int& max_time, std::vector<int>& timepoints) {
  std::ifstream time_file;
    time_file.open("Data/Time_points.dat");
    int t;
    max_time=0;
    min_time=1000;
    for (int i=0;i<1000000;i++) {
        if (!(time_file >> t)) break;
        timepoints.push_back(t);
        if (t<min_time) {
            min_time=t;
        }
        if (t>max_time) {
            max_time=t;
        }
    }
}


void OutputSummaryData (std::vector<detect>& sim_data) {
    //Matches the detections that were read in
    std::cout << "Sim data\n";
    for (int i=0;i<sim_data.size();i++) {
        std::cout << sim_data[i].day << " " << sim_data[i].cases << "\n";
    }
}

void OutputPopulationDetails (run_params& p, std::vector<int>& pop_size, std::vector<int>& pop_sum, std::vector<int> t_detects, outbreak& o) {
    std::cout << "Population size\n";
    std::cout << pop_size.size() << "\n";
    for (int i=0;i<pop_size.size();i++) {
        std::cout << pop_size[i] << " ";
    }
    std::cout << "\n";

    std::cout << "Population sum\n";
    std::cout << pop_sum.size() << "\n";
    for (int i=0;i<pop_sum.size();i++) {
        std::cout << pop_sum[i] << " ";
    }
    std::cout << "\n";

    std::cout << "Number of cases " << o.individuals.size() << "\n";
    std::cout << "Number of detections " << t_detects.size() << "\n";
    std::cout << "First detection " << o.time_first_detect << " of case symptomatic at time " << o.time_first_detect-p.time_symptom_onset_to_detect << "\n";
    std::cout << "T detects\n";
    for (int i=0;i<t_detects.size();i++) {
        std::cout << t_detects[i] << " ";
    }
    std::cout << "\n";

}

void OutputRawData (int& r0val, std::vector<int>& timepoints, std::vector< std::vector<output> >& results) {
    for (int i=0;i<timepoints.size();i++) {
        std::cout << "Timepoint " << timepoints[i] << "\n";
        std::cout << "Total accepted " << results[i][r0val].accepted << "\n";
        std::cout << "Origin times\n";
        for (int j=0;j<results[i][r0val].origin_time.size();j++) {
            std::cout << results[i][r0val].origin_time[j] << " ";
        }
        std::cout << "\n";
        std::cout << "Population sizes\n";
        for (int j=0;j<results[i][r0val].current_size.size();j++) {
            std::cout << results[i][r0val].current_size[j] << " ";
        }
        std::cout << "\n";
    }
}

void OutputAcceptanceRates (std::vector<int>& timepoints, std::vector< std::vector<output> >& results) {
    for (int r0val=1;r0val<=40;r0val++) {
      std::ofstream acc_file;
      std::ostringstream convert;
        convert << r0val;
	std::string temp=convert.str();
	std::string name = "Acceptance_rate"+temp+".dat";
        acc_file.open(name.c_str());
        for (int i=0;i<timepoints.size();i++) {
            double acc=(results[i][r0val].accepted+0.)/(results[i][r0val].tested+0.);
            acc_file << timepoints[i] << " " << acc << "\n";
        }
    }
}

void OutputOriginTimes (std::vector<int>& timepoints, std::vector< std::vector<output> >& results) {
    for (int r0val=1;r0val<=40;r0val++) {
      std::ofstream init_file;
      std::ostringstream convert;
        convert << r0val;
	std::string temp=convert.str();
	std::string name = "Origin_times"+temp+".dat";
        init_file.open(name.c_str());
        for (int i=0;i<timepoints.size();i++) {
            init_file << timepoints[i] << " ";
            for (int j=0;j<results[i][r0val].origin_time.size();j++) {
                init_file << results[i][r0val].origin_time[j] << " ";
            }
            init_file << "\n";
        }
        init_file.close();
    }
}

void OutputPopulationSizes (std::vector<int>& timepoints, std::vector< std::vector<output> >& results) {
    for (int r0val=1;r0val<=40;r0val++) {
      std::ofstream size_file;
      std::ostringstream convert;
        convert << r0val;
	std::string temp=convert.str();
	std::string name = "Current_size"+temp+".dat";
        size_file.open(name.c_str());
        for (int i=0;i<timepoints.size();i++) {
            size_file << timepoints[i] << " ";
            for (int j=0;j<results[i][r0val].current_size.size();j++) {
                if (results[i][r0val].current_size[j]>0) {
                    size_file << results[i][r0val].current_size[j] << " ";
                }
            }
            size_file << "\n";
        }
        size_file.close();
    }
}

void OutputProbabilityEnded (std::vector<int>& timepoints, std::vector< std::vector<output> >& results) {
    for (int r0val=1;r0val<=40;r0val++) {
      std::ofstream dead_file;
      std::ostringstream convert;
        convert << r0val;
	std::string temp=convert.str();
	std::string name = "End_prob"+temp+".dat";
        dead_file.open(name.c_str());
        for (int i=0;i<timepoints.size();i++) {
            double acc=(results[i][r0val].dead+0.)/(results[i][r0val].accepted+0.);
            dead_file << timepoints[i] << " " << acc << "\n";
        }
        dead_file.close();
    }
}

void OutputOutbreakDeathStatistics (run_params& p, std::vector<int>& timepoints, std::vector< std::vector<output> >& results) {
  std::ofstream pd_file;
    pd_file.open("P_Outbreak_End.dat");
    for (int i=0;i<timepoints.size();i++) { //For each of the time points
        //Make acceptance std::vector for this timepoint
        std::vector<double> acceptance;
        CalculateAcceptance (p,i,results,acceptance);
        std::vector<double> pdead;
        for (p.r0=0.1;p.r0<=4.01;p.r0=p.r0+0.1) {
            int r0val=floor((p.r0+0.001)*10);
            double pd=(results[i][r0val].dead+0.)/(results[i][r0val].accepted+0.);
            pdead.push_back(pd);
        }
        double prob=0;
        for (int j=0;j<acceptance.size();j++) {
            prob=prob+acceptance[j]*pdead[j];
        }
        pd_file << timepoints[i] << " " << prob << "\n";
    }
    pd_file.close();
}

void OutputOutbreakTimeStatistics (run_params& p, std::vector<int>& timepoints, std::vector< std::vector<output> >& results) {
    for (int i=0;i<timepoints.size();i++) { //For each of the time points
        //std::cout << "Time point " << timepoints[i] << "\n";
        std::vector<double> acceptance;
        CalculateAcceptance (p,i,results,acceptance);
        //Make a count of the origin dates for each R0, weighted by the acceptance rate for that R0
        std::vector< std::vector<double> > all_origins;
        for (p.r0=0.1;p.r0<=4.01;p.r0=p.r0+0.1) {
            std::vector<double> origins;
            int r0val=floor((p.r0+0.001)*10);
            //std::cout << "R0 " << r0val << " Size " << results[i][r0val].origin_time.size() << "\n";
            for (int k=0;k<results[i][r0val].origin_time.size();k++) {
                while (-results[i][r0val].origin_time[k]>=origins.size()) {
                    origins.push_back(0);
                }
                origins[-results[i][r0val].origin_time[k]]++;
            }
            for (int j=0;j<origins.size();j++) {
                origins[j]=origins[j]*acceptance[r0val-1];
            }
            all_origins.push_back(origins);
        }
        //Compile all of these
        std::vector<double> combined_origins;
        for (int j=0;j<all_origins.size();j++) {
            for (int k=0;k<all_origins[j].size();k++) {
                while (k>=combined_origins.size()) {
                    combined_origins.push_back(0);
                }
                combined_origins[k]=combined_origins[k]+all_origins[j][k];
            }
        }
        
        double tot=0;
        for (int j=0;j<combined_origins.size();j++) {
            tot=tot+combined_origins[j];
        }
        for (int j=0;j<combined_origins.size();j++) {
            combined_origins[j]=combined_origins[j]/tot;
        }
	std::ofstream orig_file;
	std::ostringstream convert;
        convert << timepoints[i];
	std::string temp=convert.str();
	std::string name="Origin_stats_"+temp+".dat";
        orig_file.open(name.c_str());
        for (int j=0;j<combined_origins.size();j++) {
            orig_file << j << " " << combined_origins[j] << "\n";
        }
        //std::cout << "\n";
        orig_file.close();
    }
}

void OutputOutbreakPopulationStatistics (run_params& p, std::vector<int>& timepoints, std::vector< std::vector<output> >& results) {
    //Size of population conditional on not being zero
    for (int i=0;i<timepoints.size();i++) { //For each of the time points
        std::vector<double> acceptance;
        CalculateAcceptance (p,i,results,acceptance);
        //Make a count of the population sizes for each R0, weighted by acceptance rate
        std::vector< std::vector<double> > all_sizes;
        for (p.r0=0.1;p.r0<=4.01;p.r0=p.r0+0.1) {
            std::vector<double> sizes;
            int r0val=floor((p.r0+0.001)*10);
            //std::cout << "R0 " << r0val << " Size " << results[i][r0val].current_size.size() << "\n";
            for (int k=0;k<results[i][r0val].current_size.size();k++) {
                while (results[i][r0val].current_size[k]>=sizes.size()) {
                    sizes.push_back(0);
                }
                if (results[i][r0val].current_size[k]>0) {
                    sizes[results[i][r0val].current_size[k]]++;
                }
            }
            for (int j=0;j<sizes.size();j++) {
                sizes[j]=sizes[j]*acceptance[r0val-1];
            }
            all_sizes.push_back(sizes);
        }
        //Compile all of these
        std::vector<double> combined_sizes;
        for (int j=0;j<all_sizes.size();j++) {
            for (int k=0;k<all_sizes[j].size();k++) {
                while (k>=combined_sizes.size()) {
                    combined_sizes.push_back(0);
                }
                combined_sizes[k]=combined_sizes[k]+all_sizes[j][k];
            }
        }
        double tot=0;
        for (int j=0;j<combined_sizes.size();j++) {
            tot=tot+combined_sizes[j];
        }
        for (int j=0;j<combined_sizes.size();j++) {
            combined_sizes[j]=combined_sizes[j]/tot;
        }
	std::ofstream size_file;
	std::ostringstream convert;
        convert << timepoints[i];
	std::string temp=convert.str();
	std::string name="Sizes_time_"+temp+".dat";
        size_file.open(name.c_str());
        for (int j=0;j<combined_sizes.size();j++) {
            size_file << j << " " << combined_sizes[j] << "\n";
        }
        size_file.close();
    }
}
