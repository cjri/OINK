
#include <gtest/gtest.h>
#include <vector>
#include "basicmodel.h"
#include "utilities.h"

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
    p.infection_length=5;
    std::vector<int> new_symptomatic(10, 0); // 10 time points, no new infections
    std::vector<int> total_active(10, 0); // Pre-initialized vector

    MakePopulationSize(p, new_symptomatic, total_active);

    for (const auto& active_cases : total_active) {
        EXPECT_EQ(active_cases, 0);
    }
}

TEST(MakePopulationSizeTest, ConstantNewInfections) {
    run_params p; // Infection lasts for 3 days
    p.infection_length=3;
    std::vector<int> new_symptomatic(5, 2); // 2 new infections at each of 5 time points
    std::vector<int> total_active(5, 0);

    MakePopulationSize(p, new_symptomatic, total_active);

    // Expected total_active should be: [2, 4, 6, 6, 6]
    std::vector<int> expected = {2, 4, 6, 6, 6};
    EXPECT_EQ(total_active, expected);
}

TEST(MakePopulationSizeTest, VariableNewInfections) {
    run_params p; // Infection lasts for 2 days
    p.infection_length=2;
    std::vector<int> new_symptomatic = {1, 3, 0, 2}; // Variable new infections
    std::vector<int> total_active(4, 0);

    MakePopulationSize(p, new_symptomatic, total_active);

    // Expected total_active should be: [1, 4, 3, 2]
    std::vector<int> expected = {1, 4, 3, 2};
    EXPECT_EQ(total_active, expected);
}

TEST(MakePopulationSizeTest, LongerInfectionPeriod) {
    run_params p; // Infection lasts for 6 days, longer than the simulation length
    p.infection_length=6;
    std::vector<int> new_symptomatic = {5, 0, 0}; // New infections only at the first time point
    std::vector<int> total_active(3, 0);

    MakePopulationSize(p, new_symptomatic, total_active);

    // Expected total_active should be: [5, 5, 5], as the infection affects all remaining days
    std::vector<int> expected = {5, 5, 5};
    EXPECT_EQ(total_active, expected);
}

TEST(MakePopulationSizeTest, SingleDayInfection) {
    run_params p; // Infection affects for 1 day only
    p.infection_length=1;
    std::vector<int> new_symptomatic = {1, 2, 3}; // New infections each day
    std::vector<int> total_active(3, 0);

    MakePopulationSize(p, new_symptomatic, total_active);

    // Expected total_active should directly mirror new_symptomatic
    EXPECT_EQ(total_active, new_symptomatic);
}


/*
output make_output(int tested, int accepted, int dead=0)
{
    output o;
    o.tested = tested;
    o.accepted = accepted;
    o.dead = dead;
}
*/

TEST(CalculateAcceptanceTest, UniformAcceptanceRates) {
    run_params p;
    p.R0_vals = {2.0, 2.5, 3.0}; // Example R0 values

    int time_point_index = 0; // Index of the time point for which acceptance rates are calculated
    // output = {tested, accepted, dead, current_size[], origin_time[]}
    std::vector<std::vector<output>> results = {
        {{100, 10, 0, {}, {}}, {200, 20, 0, {}, {}}, {300, 30, 0, {}, {}}} // Example output for 3 R0 values at a single time point
    };

    std::vector<double> acceptance;

    CalculateAcceptance(p, time_point_index, results, acceptance);

    // Since all R0 values have an acceptance rate of 0.1, after normalization, each should have equal weight
    ASSERT_EQ(acceptance.size(), 3);
    for (const auto& rate : acceptance) {
        EXPECT_NEAR(rate, 1.0 / 3, 0.0001);
    }
}

TEST(CalculateAcceptanceTest, DiverseAcceptanceRates) {
    run_params p;
    p.R0_vals = {2.0, 2.5, 3.0}; // Example R0 values

    int time_point_index = 0; // Index of the time point for which acceptance rates are calculated
    std::vector<std::vector<output>> results = {
       {{100, 10, 0, {}, {}}, {100, 50, 0, {}, {}}, {100, 40, 0, {}, {}}} // Example output with different acceptance rates
    };

    std::vector<double> acceptance;

    CalculateAcceptance(p, time_point_index, results, acceptance);

    // Calculate expected values for manual verification
    double expectedTotal = (10.0 / 100) + (50.0 / 100) + (40.0 / 100);
    std::vector<double> expectedRates = {
        (10.0 / 100) / expectedTotal,
        (50.0 / 100) / expectedTotal,
        (40.0 / 100) / expectedTotal
    };

    ASSERT_EQ(acceptance.size(), 3);
    for (size_t i = 0; i < acceptance.size(); ++i) {
        EXPECT_NEAR(acceptance[i], expectedRates[i], 0.0001);
    }
}

TEST(MakeRelativeTimeTest, BasicFunctionalitySortedInput) {
    std::vector<int> t_detects = {2, 4, 6, 9}; // Absolute detection times
    std::vector<int> t_detects_relative; // Vector to store relative detection times

    MakeRelativeTime(t_detects, t_detects_relative);

    std::vector<int> expected_relative_times = {0, 2, 4, 7}; // Expected relative detection times
    EXPECT_EQ(t_detects_relative, expected_relative_times);
}


TEST(CheckTerminationTest, TerminationDueToExcessiveZeroSymptoms) {
    run_params p{.infection_length = 5, .verb = 0}; 
    int t = 10; // Current simulation time
    int zeros = 7; // Consecutive zeros exceeding infection length + 1 (so infection_length+1 days with no new symptomatic cases)
    int first_evaluation_time = 3;
    int last_evaluation_time = 20;
    unsigned long n_detections = 5;
    std::vector<int> t_detects = {1, 2, 4, 6, 8}; // Arbitrary detection times
    outbreak o{.time_first_detect = 1}; // Example outbreak with first detection at time 1

    int result = CheckTermination(p, t, zeros, first_evaluation_time, last_evaluation_time, n_detections, t_detects, o);
    EXPECT_EQ(result, 1); // Expect termination
}

TEST(CheckTerminationTest, TerminationAtLastEvaluationTime) {
    run_params p{.infection_length = 5, .verb = 0}; // Example parameters
    int t = 21; // Current simulation time at last evaluation time past first detection
    int zeros = 2; // Arbitrary number of zeros not triggering termination
    int first_evaluation_time = 3;
    int last_evaluation_time = 20; // Last evaluation time
    unsigned long n_detections = 2;
    std::vector<int> t_detects = {1, 5}; // Detection times
    outbreak o{.time_first_detect = 1}; // First detection at time 1

    int result = CheckTermination(p, t, zeros, first_evaluation_time, last_evaluation_time, n_detections, t_detects, o);
    EXPECT_EQ(result, 1); // Expect termination
}

TEST(CheckTerminationTest, ContinuationWithinEvaluationRange) {
    run_params p{.infection_length = 5, .verb = 0}; // Example parameters
    int t = 11; // Current simulation time within evaluation range
    int zeros = 3; // Number of zeros not exceeding infection length + 1
    int first_evaluation_time = 5;
    int last_evaluation_time = 25; // Last evaluation time after current time
    unsigned long n_detections = 5; // Total number of detections
    std::vector<int> t_detects = {2, 4, 6, 8, 10, 12}; // Detections within the evaluation time range
    outbreak o{.time_first_detect = 2}; // Example outbreak with first detection at time 2

    int result = CheckTermination(p, t, zeros, first_evaluation_time, last_evaluation_time, n_detections, t_detects, o);
    EXPECT_EQ(result, 0); // Expect no termination
}

TEST(CheckTerminationTest, TerminationDueToExcessiveDetectionsBeforeTime) {
    run_params p{.infection_length = 5, .verb = 0}; // Example parameters
    int t = 2; // Current simulation time
    int zeros = 2; // Arbitrary number of zeros
    int first_evaluation_time = 3;
    int last_evaluation_time = 20; // Last evaluation time before current time
    unsigned long n_detections = 3; // Expected total number of detections
    std::vector<int> t_detects = {1, 2, 3, 4, 5}; // Excessive detections before the first timepoint
    outbreak o{.time_first_detect = 1}; // Example outbreak with first detection at time 1

    int result = CheckTermination(p, t, zeros, first_evaluation_time, last_evaluation_time, n_detections, t_detects, o);
    EXPECT_EQ(result, 1); // Expect termination due to excessive early detections
}

TEST(EvaluateOutbreakTest, PerfectMatchBetweenSimulationAndDetection) {
    run_params p{.infection_length = 5, .verb = 0}; // Assume this has necessary fields set
    unsigned long r0val = 0; // Index of R0 being evaluated
    std::vector<int> t_detects_relative = {0, 2, 2, 3}; // Relative detection times from the simulation
    std::vector<int> timepoints = {1, 2, 4}; // Time points at which the outbreak is evaluated
    std::vector<int> total_active_infected = {10, 20, 30}; // Number of active infections at each time point
    std::vector<detect> detections = {{0, 1}, {2, 2}, {3, 1}}; // Actual detection data
    outbreak o{.time_first_detect=2}; // Assume this has necessary fields set
    std::vector<std::vector<output>> results(timepoints.size(), std::vector<output>(1)); // Prepare results vector
    //std::vector<int> expected_accepted = {1, 1, 1}; // Expected relative detection times


    EvaluateOutbreak(p, r0val, t_detects_relative, timepoints, total_active_infected, detections, o, results);

    for (size_t i = 0; i < timepoints.size(); ++i) {
        EXPECT_EQ(results[i][r0val].tested, 1); // Each time point should be tested once
        EXPECT_EQ(results[i][r0val].accepted, 1); // Each time point should be accepted
    }
}

TEST(EvaluateOutbreakTest, NoMatchBetweenSimulationAndDetection) {
    run_params p{.infection_length = 5, .verb = 0}; 
    unsigned long r0val = 0; // Index of R0 being evaluated
    std::vector<int> t_detects_relative = {0, 0, 2, 4}; // Relative detection times from the simulation, not matching detections
    std::vector<int> timepoints = {1, 2, 3}; // Time points at which the outbreak is evaluated
    std::vector<int> total_active_infected = {0, 0, 0}; // Number of active infections at each time point, not matching detections
    std::vector<detect> detections = {{0, 1}, {2, 2}, {3, 1}}; // Actual detection data
    outbreak o{.time_first_detect=1};
    std::vector<std::vector<output>> results(timepoints.size(), std::vector<output>(1)); // Prepare results vector

    EvaluateOutbreak(p, r0val, t_detects_relative, timepoints, total_active_infected, detections, o, results);

    for (size_t i = 0; i < timepoints.size(); ++i) {
        EXPECT_EQ(results[i][r0val].tested, 1); // Each time point should be tested once
        EXPECT_EQ(results[i][r0val].accepted, 0); // No time points should be accepted
    }
}

TEST(EvaluateOutbreakTest, PartialMatchBetweenSimulationAndDetection) {
    run_params p{.infection_length = 5, .verb = 0}; 
    unsigned long r0val = 0; // Index of R0 being evaluated
    std::vector<int> t_detects_relative = {0, 2, 2, 3}; // Relative detection times from the simulation
    std::vector<int> timepoints = {1, 2, 4}; // Time points at which the outbreak is evaluated
    std::vector<int> total_active_infected = {10, 20, 30}; // Number of active infections at each time point
    std::vector<detect> detections = {{0, 1}, {2, 2}, {3,2}}; // Actual detection data
    outbreak o{.time_first_detect=1}; // Assume this has necessary fields set
    std::vector<std::vector<output>> results(timepoints.size(), std::vector<output>(1)); // Prepare results vector
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
