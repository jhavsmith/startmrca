#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
/* R CMD SHLIB cprob.c */
/* R CMD SHLIB cprobback.c */

/* paul@stat.washington.edu */
/* code  by matthew stephens */


#include <stdio.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(init_key, key_length)
unsigned long init_key[], key_length;
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 


/*
// integer generator, from 0 to (m-1) according to prob measure p
// alias method (Ripley, 1987, P232)
//
// set up the alias table
*/
void setalias(double *p, int m, int * a,double * q)
{
  int i,s;
	static int w[101];
	int nn=0;
	int np=m+1;
	for(i=1;i<(m+1);i++)
	{
		q[i]=m*p[i-1];
		if(q[i]<0.99999999999) 
		{
			w[++nn]=i;
		}
		else
		{
			w[--np]=i;
		}
	}
	for(s=1;s<m;s++)
	{
		int i0=w[s];
		int j=w[np];
		a[i0]=j;
		q[j]+=q[i0]-1;
		if(q[j]<0.99999999999) np++;
	}
	a[w[m]]=w[m];
	for(i=1;i<(m+1);i++)
	{
		q[i]+=i-1;
	}

}
//
// generate the random integer
//
int rint2(double * p, int m)
{
	int i;
	double psum = 0;
    for(i=0; i<m; i=i+1){
	  psum =psum + p[i];
	}
	for(i = 0; i<m;i=i+1){
	  p[i] = p[i]/psum;
	}
	if(m==0)
		return 0;
	else
	{
	static double q[101];
	static int a[101];
	int result;
	setalias(p,m,a,q);

//
//	produce random number
//


	double u=m*genrand_real3();
	int i=1+(int) floor(u);
	if(u<q[i])
		result=i-1;
	else
		result=a[i]-1;
	return result;
	}
}
//


/* newh is vector of length Nloci (0 and 1)
 existing_h is matrix Nhaps by Nloci
 hapcounts is vector of length Nhaps
 lambda is vector of length Nloci-1
 pos is vector of length Nloci
 nchr is total number of chromosomes in sample (used to compute theta)
 existing_h_vec is length-(Nloci * Nhaps) array of the existing haplotype (types) */


void cprobback(int * newh, int * existing_h_vec, int *p_Nloci, int *p_Nhaps, int *
	   hapcounts, int *p_nchr,  double *p_rhobar, double * lambda,
	   double * pos, double *p_prob /* p_prob is the result */,double * log_p_prob, int *copiedhap,int * seed)
{
  
  int i, locus, previouslocus, hcopy;
  double sum, Theta, TProb;
  double ** Alpha = malloc(*p_Nloci * sizeof(double *));
  int ** existing_h = malloc(*p_Nhaps * sizeof(int *));
  int j;
  double * AlphaSum = malloc(*p_Nloci * sizeof(double));
  double * TransProb = malloc( ((*p_Nloci)-1) * sizeof(double));
  int totalcounts;
  double MutProb;
  double * AlphaPointer;
  double * OldAlphaPointer;
  double lognormconst;// for normalizing
  lognormconst=0;

  init_genrand(*seed);
  
  for(i=0 ; i < *p_Nhaps ; i++)
    existing_h[i] = malloc(*p_Nloci * sizeof(int));
  
  //* fill existing_h with existing_h_vec which was transformed to a vector
  //  column by column -- so fill this way 
  for(j = 0 ; j < *p_Nloci ; j++){
    for(i = 0 ; i < *p_Nhaps ; i++){
      existing_h[i][j] = existing_h_vec[(j * (*p_Nhaps) + i)];
    }
  }
  
  for(i=0 ; i< *p_Nloci ; i++)
    Alpha[i] = malloc(*p_Nhaps * sizeof(double));
    //* TransProb[i] is probability of copying mechanism "jumping" between
  //   loci i and i+1 
  
  for(locus = 0; locus < *p_Nloci - 1; locus++)
    TransProb[locus] = 1 - exp(- (pos[locus+1]-pos[locus]) * (*p_rhobar) *
				lambda[locus]/ (*p_nchr) );

				//* Theta is given in Li and Stephens 

  for(i = 1; i < *p_nchr; i++){
    sum = sum + 1.0/i;
  }
  Theta = 1/sum;

  //* compute sum of hapcounts 
  totalcounts = 0;
  for(i=0; i< (*p_Nhaps) ; i++)
    totalcounts = totalcounts + hapcounts[i];
  MutProb = Theta/(totalcounts + Theta);
	 
  AlphaPointer = Alpha[0];
  AlphaSum[0] = 0;

  for(hcopy = 0; hcopy < *p_Nhaps; hcopy++){	
    (*AlphaPointer) = hapcounts[hcopy] * 1.0/totalcounts;
    (*AlphaPointer) *= (1-MutProb) * (newh[0] == existing_h[hcopy][0]) + MutProb * 0.5;
    AlphaSum[0] += (*AlphaPointer);
    AlphaPointer++;
  }
    lognormconst += log(AlphaSum[0]);
    log_p_prob[0] = lognormconst;

    for(hcopy = 0; hcopy < *p_Nhaps; hcopy = hcopy + 1){    
		Alpha[0][hcopy]/=AlphaSum[0];
	}    
	AlphaSum[0]=1;

  previouslocus = 0;
  //* now compute each of the other alphas in forwards algorithm
  for(locus = 1; locus < *p_Nloci; locus++){ 
    AlphaSum[locus] = 0; //* stores sum of Alpha[locus] 
    TProb = TransProb[previouslocus]; //* prob of transition 
    AlphaPointer = Alpha[locus];
    OldAlphaPointer = Alpha[previouslocus];

    for(hcopy = 0; hcopy < *p_Nhaps; hcopy = hcopy + 1){	
      *AlphaPointer = ((1-TProb) * (*OldAlphaPointer)  + TProb * AlphaSum[previouslocus] * hapcounts[hcopy] * 1.0 / totalcounts);
      *AlphaPointer *= (1-MutProb) * (newh[locus] == existing_h[hcopy][locus]) + MutProb * 0.5;
	   AlphaSum[locus] = AlphaSum[locus] + (*AlphaPointer);
	   OldAlphaPointer ++;
	   AlphaPointer ++;
    }
    lognormconst += log(AlphaSum[locus]);
    log_p_prob[locus] = lognormconst;

    for(hcopy = 0; hcopy < *p_Nhaps; hcopy = hcopy + 1){    
		Alpha[locus][hcopy]/=AlphaSum[locus];
	}    
	AlphaSum[locus]=1;
    previouslocus = locus;
  } 
	
//this is the conditional prob 
  *p_prob = lognormconst; 

//now the backwards simulation algorithm
  
  for(locus = *p_Nloci - 1; locus>0; locus--){    
    copiedhap[locus] = rint2(Alpha[locus],*p_Nhaps);
    AlphaPointer = Alpha[locus-1];
    TProb = TransProb[locus-1]; // prob of transition
    AlphaSum[locus-1] = 0;

    int ptr = 0;
    for(hcopy = 0; hcopy < *p_Nhaps; hcopy = hcopy + 1){	
		double transprob = TProb * hapcounts[hcopy] * 1.0/totalcounts;
		if(ptr == copiedhap[locus])
			transprob += (1-TProb);
		*AlphaPointer = (*AlphaPointer/AlphaSum[locus]) * transprob;
		AlphaSum[locus-1] += *AlphaPointer;
		AlphaPointer++;
		ptr ++;
	}
  }


  copiedhap[0] = rint2(Alpha[0],*p_Nhaps);
  
  //* want to add some memory free-ing commands?
  for(i=0 ; i< *p_Nloci ; i++)
    free(Alpha[i]);
  for(i=0 ; i < *p_Nhaps ; i++)
    free(existing_h[i]);
  
 
  free(Alpha);
  free(AlphaSum);
  free(TransProb);
  free(existing_h);

}


