#include "basicmodel.h"
#include "io.h"
#include "pseudorandom.h"
#include <string>


void GetRoots (const int& N, vector<int>& prs, vector< vector<int> >& all_roots) {
    GetPrimesList (N,prs);
    int tot_perms=0; //Total permutations generated
    cout << "Number of primes " << prs.size() << "\n";
    for (int i=0;i<prs.size();i++) {
        int pr=prs[i];
        int phi=prs[i]-1;
        vector<int> factors;
        GetPrimeFactors(phi,factors);
        vector<int> roots;
        GetPrimitiveRoots(N,pr,phi,factors,roots);
        all_roots.push_back(roots);
        tot_perms=tot_perms+roots.size();
    }
}

void GetPrimesList (int N_b, vector<int>& prs) {
    int lp=0;
    int val=N_b+1;
    int sqv=floor(sqrt(val+0.1))+1;
    while (lp<10) { //N.B. In practice, currently only using the first value.  Could extend this.
        int isp=1;
        for (int i=2;i<=sqv;i++) {
            if (val%i==0) {
                isp=0;
                break;
            }
        }
        if (isp==1) {
            lp++;
            prs.push_back(val);
        }
        val++;
    }
}

void GetPrimeFactors (int phi, vector<int>& factors) {
    int found=0;
    while (phi%2==0) {
        phi=phi/2;
        if (found==0) {
            factors.push_back(2);
            found=1;
        }
    }
    int sqv=floor(sqrt(phi+0.1))+1;
    for (int i=3;i<=sqv;i=i+2) {
        found=0;
        while (phi%i==0) {
            phi=phi/i;
            if (found==0) {
                factors.push_back(i);
                found=1;
            }
        }
    }
    if (phi>1) {
        factors.push_back(phi);
    }
}

int GetPrimitiveRoot (int pr, int phi, vector<int> factors) {
    for (int i=2;i<=phi;i++) {
        int found=1;
        for (int j=0;j<factors.size();j++) {
            int r=FindPower(i,phi/factors[j],pr);
            if (r==1) {
                found=0;
                break;
            }
        }
        if (found==1) {
            return i;
        }
    }
    return -1;
}

void GetPrimitiveRoots (int N, int pr, int phi, vector<int> factors, vector<int>& roots) {
    //Here set minumum root size to 1000: Prevents early-stage effects.
    for (int i=1000;i<=N-1000;i++) {
        int found=1;
        for (int j=0;j<factors.size();j++) {
            long long ii=i;
            long long jj=phi/factors[j];
            long long pp=pr;
            int r=FindPowerLL(ii,jj,pp);
            if (r==1) {
                found=0;
                break;
            }
        }
        if (found==1) {
            roots.push_back(i);
        }
    }
}

int FindPower (int x, int y, int p) {
    int res=1;
    x=x%p;
    while (y>0) {
        if (y & 1)
            res=(res*x)%p;
        y = y >> 1;
        x=(x*x)%p;
    }
    return res;
}

int FindPowerLL (long long x, long long y, long long p) {
    int res=1;
    x=x%p;
    while (y>0) {
        if (y & 1)
            res=(res*x)%p;
        y = y >> 1;
        x=(x*x)%p;
    }
    return res;
}



void SetupRandomParameters(run_params& p, vector<int>& new_infect, vector<int>& new_incubate, vector<int>& new_detect) {
    //Set up random number generation for discrete generation
    vector<double> m;
    vector<int> infect;
    vector<int> incubate;
    int N=1000000;
    
    if (p.species=="Flu") {
        ParametersFlu (p,N,m,infect,incubate);
    }
    
    //Time to infection
    int index=0;
    for (int i=0;i<infect.size();i++) {
        for (int j=0;j<infect[i];j++) {
            new_infect[index]=floor(m[i]-0.4);
            index++;
        }
    }
    //Incubation time
    index=0;
    for (int i=0;i<infect.size();i++) {
        for (int j=0;j<incubate[i];j++) {
            new_incubate[index]=floor(m[i]-0.4);
            index++;
        }
    }
    //Detection probability
    //Set up vector: Probability of detection
    //N.B. This might be time-dependent?
    for (int i=0;i<N*p.detect;i++) {
        new_detect[i]=1;
    }
}

void SetupRandomPoissonParameter (run_params& p, vector<int>& new_number) {
    int N=1000000;
    int index=0;
    for (int i=0;i<100000;i++) {
        int nvals=floor((gsl_ran_poisson_pdf(i,p.r0)*N)+0.5);
        //cout << i << " " << nvals << " " << index << "\n";
        for (int j=0;j<nvals;j++) {
            if (index<N) {
                new_number[index]=i;
                index++;
            }
        }
        if (i>=20) {
            new_number[index]=i;
            index++;
        }
        if (index==N) {
            break;
        }
    }
}

void ParametersFlu (run_params& p, int N, vector<double>& m, vector<int>& infect, vector<int>& incubate) {
    double x=0.5;
    int last_a=0;
    int last_b=0;
    int new_a=0;
    int new_b=0;
    while (new_a<N||new_b<N) {
        m.push_back(x);
        new_a=floor(gsl_cdf_weibull_P(x,p.infection_b,p.infection_a)*N+0.5);
        new_b=floor(gsl_cdf_weibull_P(x,p.incubation_b,p.incubation_a)*N+0.5);
        //cout << "Here " << x << " " << new_a << " " << new_b << "\n";
        infect.push_back(new_a-last_a);
        incubate.push_back(new_b-last_b);
        last_a=new_a;
        last_b=new_b;
        x=x+1;
    }
}

void NewPermutation (int& pp, long long& r, long long& r_orig, const vector<int>& prs, const vector< vector<int> >& all_roots, gsl_rng *rgen) {
    int rp=floor(gsl_rng_uniform(rgen)*all_roots.size());
    int rs=floor(gsl_rng_uniform(rgen)*all_roots[rp].size());
    pp=prs[rp]; //Prime number
    r=all_roots[rp][rs]; //Primitive root
    r_orig=r;
}

int GetRandomDigit (const long long& r_orig, const int& pp, const int& N, long long& r) {
    int found=0;
    while (found==0) {
        r=(r*r_orig)%pp;
        if (r<0) {
            r=r+pp;
        }
        if (r<=N) {  //Some numbers generated in prime permutation will be too large.  Remove these
            found=1;
        }
    }
    int s=r-1;
    return s;
}
