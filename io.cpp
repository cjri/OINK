#include "basicmodel.h"
#include "io.h"
#include "utilities.h"
#include <string>

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
    p.r0=1;
    p.incubation_a=7.402580682098853;
    p.incubation_b=1.7375081458523993;
    p.infection_a=1.03138989929643;
    p.infection_b=1.0025194828961315;
    p.max_infections=500000;
    p.add_limit=1000;
    p.detect=0.1;
    p.first_detect=-1; //Probability of first detection.  By default set to p.detect later.
    p.symptom_to_detect=18;
    p.replicas=1000000;
    p.run_fast=0;
    p.resolution=1; //Measured in days.  Could be changed to weekly resolution?
    p.infection_length=7;
    p.species="Flu";
    p.max_R0=4.0;
    p.test=0;
    p.more_stats=0;
    p.verb=0;
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
            p.detect=atof(argv[x]);
        } else if (p_switch.compare("--first_detect")==0) {
            x++;
            p.first_detect=atof(argv[x]);
        } else if (p_switch.compare("--resolution")==0) {
            x++;
            p.resolution=atoi(argv[x]);
        } else if (p_switch.compare("--run_fast")==0) {
            x++;
            p.run_fast=atoi(argv[x]);
        } else if (p_switch.compare("--test")==0) {
            x++;
            p.test=atoi(argv[x]);
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

        } else {
			cout << "Incorrect usage " << argv[x] << "\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
    if (p.first_detect==-1) {
        p.first_detect=p.detect;
    }
}

void ImportDetections(int& n_detections, vector<detect>& detections) {
    ifstream detect_file;
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

void ImportTimePoints(int& min_time, int& max_time, vector<int>& timepoints) {
    ifstream time_file;
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


/*void OutputCaseData (int& first_detect, ofstream case_file, vector<pat>& pdat) {
    vector<int> cases;
    for (int i=0;i<14;i++) {
        cases.push_back(0);
    }
    for (int i=0;i<pdat.size();i++) {
        int rel_time=pdat[i].time_i-first_detect; //Time of infection relative to first detection
        for (int j=0;j<14;j++) {
            if (j<=rel_time&&j>rel_time-1) { //j is first day of infection
                for (int k=j;k<j+8;k++) {
                    if (k>=0&&k<14) {
                        cases[k]++;
                    }
                }
            }
        }
    }
    case_file << "Cases ";
    for (int i=0;i<cases.size();i++) {
        case_file << cases[i] << " ";
    }
    case_file << "\n";

}*/

void OutputSummaryData (vector<detect>& sim_data) {
    //Matches the detections that were read in
    cout << "Sim data\n";
    for (int i=0;i<sim_data.size();i++) {
        cout << sim_data[i].day << " " << sim_data[i].cases << "\n";
    }
}

void OutputPopulationDetails (run_params& p, vector<int>& pop_size, vector<int>& pop_sum, vector<int> t_detects, outbreak& o) {
    cout << "Population size\n";
    cout << pop_size.size() << "\n";
    for (int i=0;i<pop_size.size();i++) {
        cout << pop_size[i] << " ";
    }
    cout << "\n";

    cout << "Population sum\n";
    cout << pop_sum.size() << "\n";
    for (int i=0;i<pop_sum.size();i++) {
        cout << pop_sum[i] << " ";
    }
    cout << "\n";

    cout << "Number of cases " << o.indiv.size() << "\n";
    cout << "Number of detections " << t_detects.size() << "\n";
    cout << "First detection " << o.first_detect << " of case symptomatic at time " << o.first_detect-p.symptom_to_detect << "\n";
    cout << "T detects\n";
    for (int i=0;i<t_detects.size();i++) {
        cout << t_detects[i] << " ";
    }
    cout << "\n";

}

void OutputRawData (int& r0val, vector<int>& timepoints, vector< vector<output> >& results) {
    for (int i=0;i<timepoints.size();i++) {
        cout << "Timepoint " << timepoints[i] << "\n";
        cout << "Total accepted " << results[i][r0val].accepted << "\n";
        cout << "Origin times\n";
        for (int j=0;j<results[i][r0val].origin_time.size();j++) {
            cout << results[i][r0val].origin_time[j] << " ";
        }
        cout << "\n";
        cout << "Population sizes\n";
        for (int j=0;j<results[i][r0val].current_size.size();j++) {
            cout << results[i][r0val].current_size[j] << " ";
        }
        cout << "\n";
    }
}

void OutputAcceptanceRates (run_params& p, vector<int>& timepoints, vector< vector<output> >& results) {
    int top=floor((p.max_R0*10)+0.5);
    for (int r0val=1;r0val<=top;r0val++) {
        ofstream acc_file;
        ostringstream convert;
        convert << r0val;
        string temp=convert.str();
        string name = "Acceptance_rate"+temp+".dat";
        acc_file.open(name.c_str());
        for (int i=0;i<timepoints.size();i++) {
            double acc=(results[i][r0val].accepted+0.)/(results[i][r0val].tested+0.);
            acc_file << timepoints[i] << " " << acc << "\n";
        }
    }
}

void OutputOriginTimes (run_params& p, vector<int>& timepoints, vector< vector<output> >& results) {
    int top=floor((p.max_R0*10)+0.5);
    for (int r0val=1;r0val<=top;r0val++) {
        ofstream init_file;
        ostringstream convert;
        convert << r0val;
        string temp=convert.str();
        string name = "Origin_times"+temp+".dat";
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

void OutputPopulationSizes (run_params& p, vector<int>& timepoints, vector< vector<output> >& results) {
    int top=floor((p.max_R0*10)+0.5);
    for (int r0val=1;r0val<=top;r0val++) {
        ofstream size_file;
        ostringstream convert;
        convert << r0val;
        string temp=convert.str();
        string name = "Current_size"+temp+".dat";
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

void OutputProbabilityEnded (run_params& p, vector<int>& timepoints, vector< vector<output> >& results) {
    int top=floor((p.max_R0*10)+0.5);
    for (int r0val=1;r0val<=top;r0val++) {
        ofstream dead_file;
        ostringstream convert;
        convert << r0val;
        string temp=convert.str();
        string name = "End_prob"+temp+".dat";
        dead_file.open(name.c_str());
        for (int i=0;i<timepoints.size();i++) {
            double acc=(results[i][r0val].dead+0.)/(results[i][r0val].accepted+0.);
            dead_file << timepoints[i] << " " << acc << "\n";
        }
        dead_file.close();
    }
}

void OutputOutbreakDeathStatistics (run_params& p, vector<int>& timepoints, vector< vector<output> >& results) {
    ofstream pd_file;
    pd_file.open("P_Outbreak_End.dat");
    for (int i=0;i<timepoints.size();i++) { //For each of the time points
        //Make acceptance vector for this timepoint
        vector<double> acceptance;
        CalculateAcceptance (p,i,results,acceptance);
        vector<double> pdead;
        for (p.r0=0.1;p.r0<=p.max_R0+0.01;p.r0=p.r0+0.1) {
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

void OutputOutbreakTimeStatistics (run_params& p, vector<int>& timepoints, vector< vector<output> >& results) {
    for (int i=0;i<timepoints.size();i++) { //For each of the time points
        //cout << "Time point " << timepoints[i] << "\n";
        vector<double> acceptance;
        CalculateAcceptance (p,i,results,acceptance);
        //Make a count of the origin dates for each R0, weighted by the acceptance rate for that R0
        vector< vector<double> > all_origins;
        for (p.r0=0.1;p.r0<=p.max_R0+0.01;p.r0=p.r0+0.1) {
            vector<double> origins;
            int r0val=floor((p.r0+0.001)*10);
            //cout << "R0 " << r0val << " Size " << results[i][r0val].origin_time.size() << "\n";
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
        vector<double> combined_origins;
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
        ofstream orig_file;
        ostringstream convert;
        convert << timepoints[i];
        string temp=convert.str();
        string name="Origin_stats_"+temp+".dat";
        orig_file.open(name.c_str());
        for (int j=0;j<combined_origins.size();j++) {
            orig_file << j << " " << combined_origins[j] << "\n";
        }
        //cout << "\n";
        orig_file.close();
    }
}

void OutputOutbreakPopulationStatistics (run_params& p, vector<int>& timepoints, vector< vector<output> >& results) {
    //Size of population conditional on not being zero
    for (int i=0;i<timepoints.size();i++) { //For each of the time points
        vector<double> acceptance;
        CalculateAcceptance (p,i,results,acceptance);
        //Make a count of the population sizes for each R0, weighted by acceptance rate
        vector< vector<double> > all_sizes;
        for (p.r0=0.1;p.r0<=p.max_R0+0.01;p.r0=p.r0+0.1) {
            vector<double> sizes;
            int r0val=floor((p.r0+0.001)*10);
            //cout << "R0 " << r0val << " Size " << results[i][r0val].current_size.size() << "\n";
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
        vector<double> combined_sizes;
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
        ofstream size_file;
        ostringstream convert;
        convert << timepoints[i];
        string temp=convert.str();
        string name="Sizes_time_"+temp+".dat";
        size_file.open(name.c_str());
        for (int j=0;j<combined_sizes.size();j++) {
            size_file << j << " " << combined_sizes[j] << "\n";
        }
        size_file.close();
    }
}
