/**
 * @file    pairwise.c
 * @brief   routine for simple pairwise epistasis test with a simple test driver
 *
 * @author  Ge Zhang <zhangge.uc@gmail.com>
 * @date    08/27/2012
 */

#include <stdio.h>
#include <math.h>

#define MATHLIB_STANDALONE
#include "Rmath.h"

double pairwise_epi_test (int cs[][3], int ct[][3]);

int main ()
{
	// boost test 0
	//int cs[3][3] = {{5,24,14}, {57,138,96}, {151,315,187}};
	//int ct[3][3] = {{4,16,10}, {71,176,86}, {120,330,200}};
	
	// boost test 1
	int cs[3][3] = {{227,225,29}, {191,189,58}, {24,33,11}};
	int ct[3][3] = {{229,219,67}, {195,199,26}, {33,42,3}};
	
	// boost test 2
	//int cs[3][3] = {{32,103,93}, {49,250,203}, {52,94,111}};
	//int ct[3][3] = {{33,121,96}, {76,203,227}, {20,124,113}};
	
	// biforce 1
	//int cs[3][3] = {{2,10,12}, {31,130,87}, {75,295,345}};
	//int ct[3][3] = {{3,14,3}, {16,90,143}, {84,353,307}};
	
	// biforce 2
	//int cs[3][3] = {{21,214,289}, {7,119,260}, {0,0,77}};
	//int ct[3][3] = {{12,123,542}, {3,69,223}, {4,14,23}};
	
	double ll, pval;
	
	// calculate log likelihood
	ll = pairwise_epi_test(cs, ct);
	// determine p-value
	pval = pchisq(ll, 4.0, 0, 0);
	
	printf("\nLog likelihood: %f; p-value: %g\n\n", ll, pval);
}


double pairwise_epi_test (int cs[][3], int ct[][3])
{
	int cn[3][3];							// two-locus genotype count in all samples
	
	int cs1[3]={0,0,0}, cs2[3]={0,0,0}; 	// genotype count of the two markers in cases
	int ct1[3]={0,0,0}, ct2[3]={0,0,0}; 	// genotype count of the two markers in controls
	int c1[3]={0,0,0}, c2[3]={0,0,0}; 		// genotype count of the first and second marker
	
	int ns=0, nt=0, n=0;					// total count in cases, controls and all sampels
	
	double pab[3][3]; 						// conditional genotype probability p(A|B)
	double pbs[3], pbt[3];					// conditional genotype probability of the second marker p(B|C)
	double psa[3], pta[3];					// conditional case/control probability given the first marker p(C|A)
	
	double ps, pt, tao=0.0;					// intermediate variables for likelihood calculation
	double ll=0.0;							// log likelihood
	
	int a, b;								// genotype of marker 1 and 2
	
	
	// Step 1. calculate marginal counts
	// NOTE: marginal counts for single marker (cs1, cs2, ct1, ct2, c1, c2, ns, nt, n) can be precalculated to increase speed
	for (a=0; a<3; a++) {
		for (b=0; b<3; b++) {
			cn[a][b] = cs[a][b] + ct[a][b];
			
			cs1[a] += cs[a][b];
			cs2[b] += cs[a][b];
			
			ct1[a] += ct[a][b];
			ct2[b] += ct[a][b];
			
			c1[a] += cn[a][b];
			c2[b] += cn[a][b];
			
		}
		
		// total count
		ns += cs1[a];
		nt += ct1[a];
		n += c1[a];
	}
	
	// Step 2. calculate (conditional) probabilities for 
	for (a=0; a<3; a++) {
		
		// p(A|B): genotype probability of the first marker conditional on the second one
		for (b=0; b<3; b++) {
			pab[a][b] = (double) cn[a][b] / c2[b];
		}
		
		// p(B|C): genotype probability of the second marker conditional on case/control 
		pbs[a] = (double) cs2[a] / ns; 
		pbt[a] = (double) ct2[a] / nt; 
		
		// p(C|A): case/control probability conditional on the genotype of first marker
		psa[a] = (double) cs1[a] / c1[a];
		pta[a] = 1.0 - psa[a]; 	// ct1[a] / c1[a];
		
	}
	
	// 3. calculate likelihood (consists of three parts)
	// NOTE: handling of cells with zero counts should be reconsidered
	for (a=0; a<3; a++) {
		for (b=0; b<3; b++) {
			// part one
			if (cs[a][b] > 0) ll += cs[a][b] * log((double) cs[a][b] / n);
			if (ct[a][b] > 0) ll += ct[a][b] * log((double) ct[a][b] / n);
			
			// part two
			ps = pab[a][b] * pbs[b] * psa[a]; 
			pt = pab[a][b] * pbt[b] * pta[a];
			
			tao += ps + pt;
			
			if (ps > 0) ll -= cs[a][b] * log(ps);
			if (pt > 0) ll -= ct[a][b] * log(pt);
		}
	}
	// part three
	ll += n*log(tao);
	
	return 2.0*ll;
}