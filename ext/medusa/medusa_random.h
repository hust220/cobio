/* 
   Generator    Relative Execution Time

   ran0     1.0
   ran1     1.3
   ran2     2.0
   ran3     0.6
   ran4     4.0

   On balance, we recommend ran1 for general use. It is portable, based on Park
   and Miller's Minimal Standard generator with an additional shuffle, and has
   no known (to us) flaws other than period exhaustion.  If you are generating
   more than 100,000,000 random numbers in a single calculation (that is, more
   than about 5% of ran1's period), we recommend the use of ran2, with its much
   longer period.

   Knuth's subtractive routine ran3 seems to be the timing winner among portable
   routines. Unfortunately the subtractive method is not so well studied, and
   not a standard. We like to keep ran3 in reserve for a "second opinion,"
   substituting it when we suspect another generator of introducing unwanted
   correlations into a calculation.

   The routine ran4 generates extremely good random deviates, and has some other
   nice properties, but it is slow. See §7.5 for discussion.  Finally, the quick
   and dirty in-line generators ranqd1 and ranqd2 are very fast, but they are
   somewhat machine dependent, and at best only as good as a 32----bit linear
   congruential generator ever is ------ in our view not good enough in many
   situations.  We would use these only in very special cases, where speed is
   critical.

   EXCEPTION FROM NUMERICAL RECIEPES
   ----------------http://www.nr.com 
   */

#pragma once

#include <cmath>
#include <iostream>
#include <cstdlib>
using namespace std;

double ran1(long& idum);
double ran2(long& idum);
double ran3(long& idum);

typedef struct IMMENSE {unsigned long l,r;} immense;
typedef struct GREAT {unsigned short l,c,r;} great;
static unsigned long bit[33];  /* defining declaration */

unsigned long getbit(immense source, int bitno, int nbits);
void ks(immense key, int n, great* kn);

void cyfun(unsigned long ir, great k,unsigned long* iout);
void des(immense inp, immense key, int* newkey, int isw, immense* out);
double ran4(long& idum);

typedef enum {_RAN1_=0, _RAN2_, _RAN3_, _RAN4_} random_t;
typedef double (*ran)(long&);

class randomGenerator{
    long idum;
    ran f;
public:
    randomGenerator() {
        init(_RAN2_, 11);
    }

    randomGenerator(random_t type, long seed = -101) {
        init(type, seed);
    }

    void init(random_t type, long seed = -101){
        idum = -labs(seed);
        switch(type){
            case _RAN1_:
                f = ran1;
                break;
            case _RAN2_:
                f = ran2;
                break;
            case _RAN3_:
                f = ran3;
                break;
            case _RAN4_:
                f = ran4;
                break;
            default:
                cerr << "Wrong input" << endl;
                exit(1);
                break;
        }
        f(idum);
    }

    ~randomGenerator(){};

    double next(){
        return f(idum);
    }

    /*------------------------------
      p(x) = Exp(-x*x/2)/sqrt(2*PI);
      -----------------------------*/
    double nextGauss(){
        static int iset=0;
        static double gset;
        double fac,r,v1,v2;

        if  (iset == 0) {
            do {
                v1=2.0*f(idum)-1.0;
                v2=2.0*f(idum)-1.0;
                r=v1*v1+v2*v2;
            } while (r >= 1.0);
            fac=sqrt(-2.0*log(r)/r);
            gset=v1*fac;
            iset=1;
            return v2*fac;
        } else {
            iset=0;
            return gset;
        }
    }

};

/*  
    int main(){
    randomGenerator r(_RAN2_, -1);
    for(int i=0; i<10000; i++)
    cout << r.nextGauss() << endl;
    }
    */
