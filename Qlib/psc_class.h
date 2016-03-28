#ifndef PSC_H
#define PSC_H

#include <cmath>
#include <string>
#include <iostream>
#include "parameters.h"
#include "functions.h"
using namespace std;


class PSC
{
    public:
        PSC(string name, double taod, double taor, double tinval=TINVAL, int totaltimestep=CLK);
        ~PSC();
        int getcurrentstep(){return currentstep;}
        int gettimeinval(){return timeinval;}
        double getg(){return g[currentstep];}
        int setcurrentstep(int s){currentstep = s;}
        void reset();
        void init(double gg, double hh){g[0]=gg;h[0]=hh;}
        bool regular_evolve();
        bool spikecome(double strong, double ts);
        ostream& operator<<(ostream &s);
        void print_state(ostream &s);
    private:
        double timeinval;
        int totalstep;
        string PSCname;
        double tao_d;
        double tao_r;
    protected:
        int currentstep;
        double *h;
        double *g;
};

PSC::PSC(string name, double taod, double taor, double tinval, int totaltimestep)
{
     PSCname = name;
     currentstep = 0;
     tao_d = taod;
     tao_r = taor;
     timeinval = tinval;
     totalstep = totaltimestep;
     h = NULL;
     g = NULL;
     if(totalstep>0)
     {
          h = new double[totalstep+1];
          g = new double[totalstep+1];
     }
}
PSC::~PSC()
{
     if(h)
         delete [] h;
     if(g)
         delete [] g;
}
bool PSC::regular_evolve()
{
    if( currentstep+1 > totalstep )
        return false;

    //h[currentstep+1] = h[currentstep]*exp(-timeinval/tao_r);
    //g[currentstep+1] = g[currentstep]*exp(-timeinval/tao_d)+h[currentstep]*tao_r/(tao_r-tao_d)*(exp(-timeinval/tao_r)-exp(-timeinval/tao_d));
    h[currentstep+1]=regular_s(h[currentstep],timeinval,tao_r);
    g[currentstep+1]=regular_g2(h[currentstep],g[currentstep],timeinval,tao_r,tao_d);
    currentstep++;
    //cout<<PSCname<<":"<<currentstep<<" evolved to "<<g[currentstep]<<endl;
    return true;
}
bool PSC::spikecome(double strong, double ts)
{
    if(currentstep<1)
        return false;
    //double t = timeinval - ts;
    //h[currentstep] += strong / tao_r * exp(-t/tao_r);
    //g[currentstep] += strong / (tao_r-tao_d) *( exp(-t/tao_r)-exp(-t/tao_d) );
    h[currentstep] = add_s(h[currentstep],timeinval,ts,strong,tao_r);
    g[currentstep] = add_g2(g[currentstep],timeinval,ts,strong,tao_r,tao_d);
    //cout<<PSCname<<":"<<currentstep<<" recieve spike"<<endl;
    return true;
}
void PSC::reset()
{
     h[0] = h[currentstep];
     g[0] = g[currentstep];
     currentstep = 0;
}
ostream& PSC::operator<<(ostream &s)
{
    for(int i=0; i<currentstep; i++)
        s<<h[currentstep+1]<<'\t';
    s<<endl;
    for(int i=0; i<currentstep; i++)
        s<<g[currentstep+1]<<'\t';
    s<<endl;
    return s;
}
void PSC::print_state(ostream &s)
{
    for(int i=0; i<currentstep; i++)
        s<<h[i+1]<<'\t';
    s<<endl;
    for(int i=0; i<currentstep; i++)
        s<<g[i+1]<<'\t';
    s<<endl;
}

#endif
