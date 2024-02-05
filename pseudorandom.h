#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

void GetRoots (const int& N, vector<int>& prs, vector< vector<int> >& all_roots);
void GetPrimesList (int N_b, vector<int>& prs);
void GetPrimeFactors (int phi, vector<int>& factors);
int GetPrimitiveRoot (int pr, int phi, vector<int> factors);
void GetPrimitiveRoots (int N, int pr, int phi, vector<int> factors, vector<int>& roots);
int FindPower (int x, int y, int p);
int FindPowerLL (long long x, long long y, long long p);

void NewPermutation (int& pp, long long& r, long long& r_orig, const vector<int>& prs, const vector< vector<int> >& all_roots, gsl_rng *rgen);
int GetRandomDigit (const long long& r_orig, const int& pp, const int& N, long long& r);


void SetupRandomParameters(run_params& p, vector<int>& new_infect, vector<int>& new_incubate, vector<int>& new_detect);
void SetupRandomPoissonParameter (run_params& p, vector<int>& new_number);

void ParametersFlu (run_params& p, int N, vector<double>& m, vector<int>& infect, vector<int>& incubate);

