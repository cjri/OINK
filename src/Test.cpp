
#include <gtest/gtest.h>
#include <vector>
#include "basicmodel.h"
#include "utilities.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>

TEST(CompareWithDataTest, DifferentAfterLastTimepoint) {
    // Assume CompareWithData is part of a namespace or adjust accordingly
    std::vector<int> timepoints = {2, 3, 5}; 
    std::vector<detect> sim_data = {{1, 100}, {2, 200}, {3, 300}, {7, 20}};
    std::vector<detect> detections = {{1, 100}, {2, 200}, {3, 300}, {6, 20}};

    // Call the CompareWithData function with the test data
    auto acceptedTimepoints = CompareWithData(timepoints, sim_data, detections);

    EXPECT_EQ(acceptedTimepoints, 3);
}

TEST(CompareWithDataTest, ExtraSimData) {
    std::vector<int> timepoints = {2, 3, 5}; 
    std::vector<detect> sim_data = {{1, 100}, {2, 200}, {3, 300}, {4, 20}};
    std::vector<detect> detections = {{1, 100}, {2, 200}, {3, 300}};

    // Call the CompareWithData function with the test data
    auto acceptedTimepoints = CompareWithData(timepoints, sim_data, detections);

    EXPECT_EQ(acceptedTimepoints, 2);
}

TEST(CompareWithDataTest, ExtraDetectionData) {
    std::vector<int> timepoints = {1, 3, 5}; 
    std::vector<detect> sim_data = {{0, 100}, {2, 200}, {3, 300}};
    std::vector<detect> detections = {{0, 100}, {2, 200}, {3, 300}, {4,20}};

    // Call the CompareWithData function with the test data
    auto acceptedTimepoints = CompareWithData(timepoints, sim_data, detections);

    EXPECT_EQ(acceptedTimepoints, 2);
}

TEST(CompareWithDataTest, ExtraDetectionDataSimAfter) {
    std::vector<int> timepoints = {1, 3, 5}; 
    std::vector<detect> sim_data = {{0, 100}, {2, 200}, {3, 300}, {7,10} };
    std::vector<detect> detections = {{0, 100}, {2, 200}, {3, 300}, {4,20}};

    // Call the CompareWithData function with the test data
    auto acceptedTimepoints = CompareWithData(timepoints, sim_data, detections);

    EXPECT_EQ(acceptedTimepoints, 2);
}

TEST(ConstructSummaryDataTest, HandlesEmptyInput) {
    std::vector<int> t_detects;
    std::vector<detect> sim_data;
    ConstructSummaryData(t_detects, sim_data);
    EXPECT_TRUE(sim_data.empty());
}

TEST(ConstructSummaryDataTest, NonEmptyInputWithoutZero) {
    std::vector<int> t_detects = {1, 1, 2};
    std::vector<detect> sim_data;
    ConstructSummaryData(t_detects, sim_data);

    ASSERT_EQ(sim_data.size(), 2);
    EXPECT_EQ(sim_data[0].day, 1);
    EXPECT_EQ(sim_data[0].cases, 2);
    EXPECT_EQ(sim_data[1].day, 2);
    EXPECT_EQ(sim_data[1].cases, 1);
}

TEST(ConstructSummaryDataTest, SingleDetection) {
    std::vector<int> t_detects = {0};
    std::vector<detect> sim_data;
    ConstructSummaryData(t_detects, sim_data);

    ASSERT_EQ(sim_data.size(), 1);
    EXPECT_EQ(sim_data[0].day, 0);
    EXPECT_EQ(sim_data[0].cases, 1);
}

TEST(ConstructSummaryDataTest, ConsecutiveDetectionDays) {
    std::vector<int> t_detects = {0, 1, 1, 2, 3, 3, 3};
    std::vector<detect> sim_data;
    ConstructSummaryData(t_detects, sim_data);

    ASSERT_EQ(sim_data.size(), 4);
    EXPECT_EQ(sim_data[0].day, 0);
    EXPECT_EQ(sim_data[0].cases, 1);
    EXPECT_EQ(sim_data[1].day, 1);
    EXPECT_EQ(sim_data[1].cases, 2);
    EXPECT_EQ(sim_data[2].day, 2);
    EXPECT_EQ(sim_data[2].cases, 1);
    EXPECT_EQ(sim_data[3].day, 3);
    EXPECT_EQ(sim_data[3].cases, 3);
}


TEST(ConstructSummaryDataTest, NonConsecutiveDetectionDays) {
    std::vector<int> t_detects = {0, 2, 2, 5};
    std::vector<detect> sim_data;
    ConstructSummaryData(t_detects, sim_data);

    ASSERT_EQ(sim_data.size(), 3);
    EXPECT_EQ(sim_data[0].day, 0);
    EXPECT_EQ(sim_data[0].cases, 1);
    EXPECT_EQ(sim_data[1].day, 2);
    EXPECT_EQ(sim_data[1].cases, 2);
    EXPECT_EQ(sim_data[2].day, 5);
    EXPECT_EQ(sim_data[2].cases, 1);
}

TEST(MakePopulationSizeTest, NoNewInfections) {
    run_params p; // Assuming infection length is 5 days
    p.infection_length_min=5;
    p.infection_length_max=5;
    std::vector<int> new_symptomatic(10, 0); // 10 time points, no new infections
    std::vector<int> total_active(10, 0); // Pre-initialized vector

    outbreak o;
    // Construct outbreak - empty
    SetupOutbreak(o);

    MakePopulationSize(p, o, total_active);

    for (const auto& active_cases : total_active) {
        EXPECT_EQ(active_cases, 0);
    }
}



TEST(MakePopulationSizeTest, ConstantNewInfections) {
    run_params p; // Infection lasts for 3 days

    std::vector<int> total_active(5, 0);

    outbreak o;
    // Construct outbreak - empty
    SetupOutbreak(o);
    for(unsigned long i =0; i<5; i++)
    {
        o.individuals.push_back( patient{.time_symptom_onset = static_cast<int>(i), .infection_length=3});
        o.individuals.push_back( patient{.time_symptom_onset = static_cast<int>(i), .infection_length=3});
    }

    MakePopulationSize(p, o, total_active);

    // Expected total_active should be: [2, 4, 6, 6, 6]
    std::vector<int> expected = {2, 4, 6, 6, 6};
    EXPECT_EQ(total_active, expected);
}
TEST(MakePopulationSizeTest, VariableNewInfections) {
    run_params p; // Infection lasts for 2 days

    unsigned long infection_length=2;

    std::vector<int> new_symptomatic = {1, 3, 0, 2}; // Variable new infections
    std::vector<int> total_active(4, 0);

    outbreak o;
    // Construct outbreak - empty
    SetupOutbreak(o);
    for(unsigned long i =0; i<new_symptomatic.size(); i++)
    {
        for(int j=0; j< new_symptomatic[i]; j++)
        o.individuals.push_back( patient{.time_symptom_onset = static_cast<int>(i), .infection_length=static_cast<int>(infection_length)});
    }

    MakePopulationSize(p, o, total_active);

    // Expected total_active should be: [1, 4, 3, 2]
    std::vector<int> expected = {1, 4, 3, 2};
    EXPECT_EQ(total_active, expected);
}

TEST(MakePopulationSizeTest, LongerInfectionPeriod) {
    run_params p; // Infection lasts for 6 days, longer than the simulation length
    unsigned long infection_length=6;
    std::vector<int> new_symptomatic = {5, 0, 0}; // New infections only at the first time point
    std::vector<int> total_active(3, 0);
    outbreak o;
    // Construct outbreak - empty
    SetupOutbreak(o);
    for(unsigned long i =0; i<new_symptomatic.size(); i++)
    {
        for(int j=0; j< new_symptomatic[i]; j++)
        o.individuals.push_back( patient{.time_symptom_onset = static_cast<int>(i), .infection_length=static_cast<int>(infection_length)});
    }

    MakePopulationSize(p, o, total_active);

    // Expected total_active should be: [5, 5, 5], as the infection affects all remaining days
    std::vector<int> expected = {5, 5, 5};
    EXPECT_EQ(total_active, expected);
}


TEST(MakePopulationSizeTest, SingleDayInfection) {
    run_params p; // Infection affects for 1 day only
    unsigned long infection_length=1;

    std::vector<int> new_symptomatic = {1, 2, 3}; // New infections each day
    std::vector<int> total_active(3, 0);
    outbreak o;
    // Construct outbreak - empty
    SetupOutbreak(o);
    for(unsigned long i =0; i<new_symptomatic.size(); i++)
    {
        for(int j=0; j< new_symptomatic[i]; j++)
        o.individuals.push_back( patient{.time_symptom_onset = static_cast<int>(i), .infection_length=static_cast<int>(infection_length)});
    }

    MakePopulationSize(p, o, total_active);

    // Expected total_active should directly mirror new_symptomatic
    EXPECT_EQ(total_active, new_symptomatic);
}

TEST(NewCaseTest, BoundDetectionTime) {
  static gsl_rng *rgen;
    run_params p;
    
    p.time_symptom_onset_to_detect_min=3;
    p.time_symptom_onset_to_detect_max=5;
    p.infection_a = 1;
    p.infection_b = 1;
    p.incubation_a = 7;
    p.incubation_b = 1.7;
    p.probability_first_detect = 1.0;
    p.probability_detect = 1.0;
    p.infection_length_min=5;
    p.infection_length_max=5;
    std::vector<int> t_detects;
    std::vector<int> t_detects_after;
    outbreak o;
    gsl_rng_env_setup();
    rgen = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rgen, 1);
    SetupOutbreak(o);
    MakeNewCase(p, -1, t_detects, t_detects_after, o, rgen);
    MakeNewCase(p, 0, t_detects, t_detects_after, o, rgen);
    EXPECT_EQ(t_detects.size(), 2);
    EXPECT_EQ(t_detects_after.size(), 0);
    int t = t_detects[0] - o.individuals[0].time_symptom_onset;
    EXPECT_LE(t, 5);
    EXPECT_GE(t, 3);
    t = t_detects[1] - o.individuals[1].time_symptom_onset;
    EXPECT_LE(t, 5);
    EXPECT_GE(t, 3);
}
TEST(NewCaseTest, NoDetection) {
  static gsl_rng *rgen;
    run_params p;
    
    p.time_symptom_onset_to_detect_min=3;
    p.time_symptom_onset_to_detect_max=5;
    p.infection_a = 1;
    p.infection_b = 1;
    p.incubation_a = 7;
    p.incubation_b = 1.7;
    p.probability_first_detect = 0.0;
    p.probability_detect = 0.0;
    p.probability_detect_enhanced = 1.0;
    p.infection_length_min=5;
    p.infection_length_max=5;
    std::vector<int> t_detects;
    std::vector<int> t_detects_after;
    outbreak o;
    gsl_rng_env_setup();
    rgen = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rgen, 1);
    SetupOutbreak(o);
    MakeNewCase(p, -1, t_detects, t_detects_after, o, rgen);
    MakeNewCase(p, 0, t_detects, t_detects_after, o, rgen);
    EXPECT_EQ(t_detects.size(), 0);
    EXPECT_EQ(t_detects_after.size(), 2);
}
    

TEST(MakeRelativeTimeTest, BasicFunctionalitySortedInput) {
    std::vector<int> t_detects = {2, 4, 6, 9}; // Absolute detection times
    std::vector<int> t_detects_relative; // Vector to store relative detection times

    MakeRelativeTime(t_detects, t_detects_relative);

    std::vector<int> expected_relative_times = {0, 2, 4, 7}; // Expected relative detection times
    EXPECT_EQ(t_detects_relative, expected_relative_times);
}
//int CheckTermination(const run_params &p, const int t, const int zeros, const int extreme_infection_time, const int first_evaluation_time, const int last_evaluation_time, const unsigned long n_detections, const std::vector<int> &t_detects, const outbreak &o);


TEST(CheckTerminationTest, TerminationDueToExcessiveZeroSymptoms) {
    run_params p{.verb = 0}; 
    int extreme_infection_time = 5;
    int t = 10; // Current simulation time
    int zeros = 7; // Consecutive zeros exceeding infection length + 1 (so infection_length+1 days with no new symptomatic cases)
    int first_evaluation_time = 3;
    int last_evaluation_time = 20;
    unsigned long n_detections = 5;
    std::vector<int> t_detects = {1, 2, 4, 6, 8}; // Arbitrary detection times
    outbreak o{.time_first_detect = 1}; // Example outbreak with first detection at time 1
    int result = CheckTermination(p, t, zeros, extreme_infection_time, first_evaluation_time, last_evaluation_time, n_detections, t_detects, o);
    EXPECT_EQ(result, 1); // Expect termination
}

TEST(CheckTerminationTest, TerminationAtLastEvaluationTime) {
    run_params p{.verb = 0}; // Example parameters
    int extreme_infection_time = 5;
    int t = 21; // Current simulation time at last evaluation time past first detection
    int zeros = 2; // Arbitrary number of zeros not triggering termination
    int first_evaluation_time = 3;
    int last_evaluation_time = 20; // Last evaluation time
    unsigned long n_detections = 2;
    std::vector<int> t_detects = {1, 5}; // Detection times
    outbreak o{.time_first_detect = 1}; // First detection at time 1
    int result = CheckTermination(p, t, zeros, extreme_infection_time, first_evaluation_time, last_evaluation_time, n_detections, t_detects, o);
    EXPECT_EQ(result, 1); // Expect termination
}

TEST(CheckTerminationTest, ContinuationWithinEvaluationRange) {
    run_params p{.infection_length_min = 5, .infection_length_max=5}; // Example parameters
    int extreme_infection_time = 5;
    int t = 11; // Current simulation time within evaluation range
    int zeros = 3; // Number of zeros not exceeding infection length + 1
    int first_evaluation_time = 5;
    int last_evaluation_time = 25; // Last evaluation time after current time
    unsigned long n_detections = 5; // Total number of detections
    std::vector<int> t_detects = {2, 4, 6, 8, 10, 12}; // Detections within the evaluation time range
    outbreak o{.time_first_detect = 2}; // Example outbreak with first detection at time 2
    int result = CheckTermination(p, t, zeros, extreme_infection_time, first_evaluation_time, last_evaluation_time, n_detections, t_detects, o);
    EXPECT_EQ(result, 0); // Expect no termination
}

TEST(CheckTerminationTest, TerminationDueToExcessiveDetectionsBeforeTime) {
    run_params p{ .verb = 0}; // Example parameters
    int extreme_infection_time = 5;
    int t = 2; // Current simulation time
    int zeros = 2; // Arbitrary number of zeros
    int first_evaluation_time = 3;
    int last_evaluation_time = 20; // Last evaluation time before current time
    unsigned long n_detections = 3; // Expected total number of detections
    std::vector<int> t_detects = {1, 2, 3, 4, 5}; // Excessive detections before the first timepoint
    outbreak o{.time_first_detect = 1}; // Example outbreak with first detection at time 1
    int result = CheckTermination(p, t, zeros, extreme_infection_time, first_evaluation_time, last_evaluation_time, n_detections, t_detects, o);
    EXPECT_EQ(result, 1); // Expect termination due to excessive early detections
}

TEST(EvaluateOutbreakTest, PerfectMatchBetweenSimulationAndDetection) {
    run_params p{.infection_length_min = 5, .infection_length_max=5, .verb = 0}; // Assume this has necessary fields set
    unsigned long r0val = 0; // Index of R0 being evaluated
    std::vector<int> t_detects_relative = {0, 2, 2, 3}; // Relative detection times from the simulation
    std::vector<int> timepoints = {1, 2, 4}; // Time points at which the outbreak is evaluated
    std::vector<int> total_active_infected = {10, 20, 30}; // Number of active infections at each time point
    std::vector<detect> detections = {{0, 1}, {2, 2}, {3, 1}}; // Actual detection data
    outbreak o{.time_first_detect=2}; // Assume this has necessary fields set
    std::vector<std::vector<output> > results(timepoints.size(), std::vector<output>(1)); // Prepare results vector
    //std::vector<int> expected_accepted = {1, 1, 1}; // Expected relative detection times


    EvaluateOutbreak(p, r0val, t_detects_relative, timepoints, total_active_infected, detections, o, results);

    for (size_t i = 0; i < timepoints.size(); ++i) {
        EXPECT_EQ(results[i][r0val].tested, 1); // Each time point should be tested once
        EXPECT_EQ(results[i][r0val].accepted, 1); // Each time point should be accepted
    }
}

TEST(EvaluateOutbreakTest, NoMatchBetweenSimulationAndDetection) {
    run_params p{.infection_length_min = 5, .infection_length_max = 5, .verb = 0}; 
    unsigned long r0val = 0; // Index of R0 being evaluated
    std::vector<int> t_detects_relative = {0, 0, 2, 4}; // Relative detection times from the simulation, not matching detections
    std::vector<int> timepoints = {1, 2, 3}; // Time points at which the outbreak is evaluated
    std::vector<int> total_active_infected = {0, 0, 0}; // Number of active infections at each time point, not matching detections
    std::vector<detect> detections = {{0, 1}, {2, 2}, {3, 1}}; // Actual detection data
    outbreak o{.time_first_detect=1};
    std::vector<std::vector<output> > results(timepoints.size(), std::vector<output>(1)); // Prepare results vector

    EvaluateOutbreak(p, r0val, t_detects_relative, timepoints, total_active_infected, detections, o, results);

    for (size_t i = 0; i < timepoints.size(); ++i) {
        EXPECT_EQ(results[i][r0val].tested, 1); // Each time point should be tested once
        EXPECT_EQ(results[i][r0val].accepted, 0); // No time points should be accepted
    }
}

TEST(EvaluateOutbreakTest, PartialMatchBetweenSimulationAndDetection) {
    run_params p{.infection_length_min = 5, .infection_length_max=5, .verb = 0}; 
    unsigned long r0val = 0; // Index of R0 being evaluated
    std::vector<int> t_detects_relative = {0, 2, 2, 3}; // Relative detection times from the simulation
    std::vector<int> timepoints = {1, 2, 4}; // Time points at which the outbreak is evaluated
    std::vector<int> total_active_infected = {10, 20, 30}; // Number of active infections at each time point
    std::vector<detect> detections = {{0, 1}, {2, 2}, {3,2}}; // Actual detection data
    outbreak o{.time_first_detect=1}; // Assume this has necessary fields set
    std::vector<std::vector<output> > results(timepoints.size(), std::vector<output>(1)); // Prepare results vector
    std::vector<int> expected_accepted = {1, 1, 0}; // Expected relative detection times

    EvaluateOutbreak(p, r0val, t_detects_relative, timepoints, total_active_infected, detections, o, results);

    for (size_t i = 0; i < timepoints.size(); ++i) {
        EXPECT_EQ(results[i][r0val].tested, 1); // Each time point should be tested once
        EXPECT_EQ(results[i][r0val].accepted, expected_accepted[i]); // Each time point should be accepted
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
