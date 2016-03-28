#ifndef NEURON_H
#define NEURON_H
#include "psc_class.h"
#include <iostream>
#include <vector>
#include <map>
#include <stdlib.h>
#include "functions.h"
#include "parameters.h"
using namespace std;

class Neuron
{
    public:
        Neuron(double phase, double theta, int simple, int excit, double Alpha=ALPHA, double tinval=TINVAL, int totaltimestep=CLK, double th = VTHR, double gL = 50);
        ~Neuron();
        void initvoltage(double v){ voltage[0] = v; }
        void initpsc(string name,double g, double h);

        void printparameters();
        void print_voltage(ostream &s);
        void print_g(ostream &s);
        void print_spiketime(ostream &s);
        void print_psc(string name, ostream &s);

        double getVoltage(){return voltage[currentstep];}
        double getGi(){return gi[currentstep];}
        double getGe(){return ge[currentstep];}

        bool evolve(double theta, double phase, double constrast);
        bool evolve();

        void Espikecome(double strong, double t);
        void Ispikecome(double strong, double t);

        int issimple(){return isSimple;}
        int isexcit(){return isExcit;}
        void reset();
        bool checkspike(double *t);
    private:
        int totalstep;
        double timeinval;//unit: second

        int isSimple;
        int isExcit;
        double pphase;
        double ptheta;
        double tao_ref;
        double delay;
        double v_th;

    protected:
        /*
         * all following variables evoluating with time
         */
        int currentstep;
        double current_time;
        /*
         * following : conductance variables
         */
        double alpha;
        map<string,PSC*> pscmap;
        PSC *nmda_psc;
        PSC *gaba_psc;
        PSC *ampa_psc;
        double glgnE;
        double glgnI;
        double *ge;
        double *gi;
        double gl;
        bool recieve_lgn_Espike(double theta, double phase, double constrast);
        bool recieve_lgn_Ispike(double fr);
        void update_g();
        /*
         * following : voltage variables
         */
        vector<double> spike_time;
        double ts;
        double *voltage;
        bool *existspike;
};
Neuron::Neuron(double phase, double theta, int simple, int excit, double Alpha, double tinval, int totaltimestep, double th, double gL)
{
     pphase = phase;
     ptheta = theta;
     isSimple = simple;
     delay = -1000;
     v_th = th;
     gl = gL;
     isExcit = excit;
     if(isExcit)
         tao_ref = REF_E;
     else
         tao_ref = REF_I;
     alpha = Alpha;
     timeinval = tinval;
     totalstep = totaltimestep;
    ge = gi = voltage = NULL;
    existspike = NULL;
    if(totalstep>0)
    {
        ge = new double[totalstep+1];
        gi = new double[totalstep+1];
        voltage = new double[totalstep+1];
        existspike = new bool[totalstep+1];
    }
    currentstep = 0;
    current_time = 0;
    ampa_psc = new PSC("ampa",AMPAD,AMPAR,timeinval,totalstep);
    nmda_psc = new PSC("nmda",NMDAD,NMDAR,timeinval,totalstep);
    gaba_psc = new PSC("gaba",GABAD,GABAR,timeinval,totalstep);
    pscmap.insert( pair<string,PSC*>("ampa",ampa_psc) );
    pscmap.insert( pair<string,PSC*>("nmda",nmda_psc) );
    pscmap.insert( pair<string,PSC*>("gaba",gaba_psc) );

}
Neuron::~Neuron()
{
    if(ge)
        delete [] ge;
    if(gi)
        delete [] gi;
    if(voltage)
        delete [] voltage;
    if(existspike)
        delete [] existspike;
    if(nmda_psc)
        delete nmda_psc;
    if(ampa_psc)
        delete ampa_psc;
    if(gaba_psc)
        delete gaba_psc;
}

void Neuron::initpsc(string name,double g, double h)
{
     map<string,PSC*>::iterator it = pscmap.find(name);
     if( it != pscmap.end() )
        it->second->init(g,h);
}
void Neuron::printparameters()
{
    cout<<voltage[0]<<endl;
}
bool Neuron::evolve()
{
    if( currentstep++ > totalstep )
        return false;
    current_time += timeinval;

    nmda_psc -> regular_evolve();
    ampa_psc -> regular_evolve();
    gaba_psc -> regular_evolve();

    update_g();

    if(delay>0) // still in refraction
    {
        voltage[currentstep] = 0;
        delay-=timeinval;
    }
    else if(delay==0) // just finish the refraction
    {
        delay = -1000;
        voltage[currentstep] = find_v2(ge[currentstep-1],ge[currentstep],gi[currentstep-1],gi[currentstep],ts,timeinval);
    }
    else //regular
    {
        //RK2;
        voltage[currentstep] = RK2(voltage[currentstep-1],ge[currentstep-1],ge[currentstep],gi[currentstep-1],gi[currentstep],timeinval);
    }
    return true;
}
bool Neuron::evolve(double theta, double phase, double constrast)
{
    if( currentstep++ > totalstep )
        return false;
    current_time += timeinval;

    nmda_psc -> regular_evolve();
    ampa_psc -> regular_evolve();
    gaba_psc -> regular_evolve();


    if(recieve_lgn_Ispike(lgnIfr))
    {
        //cout<<currentstep<<'\t'<<"Ispikecome"<<endl;
        Ispikecome(lgnI_str,timeinval/2);
    }
    if(isSimple&&recieve_lgn_Espike(theta,phase,constrast))
    {
        Espikecome(lgnE_str,timeinval/2);
        //cout<<currentstep<<'\t'<<"Espikecome"<<endl;
    }

    update_g();
    
    if(delay>0) // still in refraction
    {
        voltage[currentstep] = 0;
        delay-=timeinval;
    }
    else if(delay==0) // just finish the refraction
    {
        delay = -1000;
        voltage[currentstep] = find_v2(ge[currentstep-1],ge[currentstep],gi[currentstep-1],gi[currentstep],ts,timeinval);
    }
    else //regular
    {
        //RK2;
        voltage[currentstep] = RK2(voltage[currentstep-1],ge[currentstep-1],ge[currentstep],gi[currentstep-1],gi[currentstep],timeinval);
    }
    return true;
}
bool Neuron::recieve_lgn_Espike(double theta, double phase, double constrast)
{
    double mysin = sin(phase-pphase);
    double mycos = cos(2*(theta-ptheta));
    double f = 80*constrast*(1+0.5*mysin*mycos);
    //cout<<f<<endl;
    int r;
    r = rand()%1001;
    //cout<<"Elgn pro:"<<r<<'\t'<<1000*f*timeinval<<endl;
    if(r < 1000*f*timeinval)
        return true;
    else
        return false;
}

bool Neuron::recieve_lgn_Ispike(double fr)
{
    int r;
    r = rand()%1001;
    //cout<<"Ilgn pro:"<<r<<'\t'<<1000*fr*timeinval<<endl;
    if(r < 1000*fr*timeinval)
        return true;
    else
        return false;
}
void Neuron::update_g()
{
    ge[currentstep] = alpha*ampa_psc->getg() + (1-alpha)*nmda_psc->getg();
    gi[currentstep] = gaba_psc->getg();
    //cout<<"update g, timestep:"<<currentstep<<" ge: "<<ge[currentstep]<<" gi: "<<gi[currentstep]<<endl;
}
void Neuron::Ispikecome(double strong, double t)
{
    gaba_psc -> spikecome(strong,t);
    //update_g();
}
void Neuron::Espikecome(double strong, double t)
{
    ampa_psc -> spikecome(strong,t);
    nmda_psc -> spikecome(strong,t);
    //update_g();
}
bool Neuron::checkspike(double *t)
{
    if(voltage[currentstep]<=v_th)
        return false;
    ts = find_ts(v_th,voltage[currentstep-1],voltage[currentstep],timeinval);
    *t = ts;
    voltage[currentstep] = 0;
    double time = ts + current_time - timeinval;
    spike_time.push_back(time);
    delay = tao_ref;
    return true;
}
void Neuron::reset()
{
    gi[0] = gi[currentstep];
    ge[0] = ge[currentstep];
    voltage[0] = voltage[currentstep];
    currentstep = 0;
    ampa_psc->reset();
    nmda_psc->reset();
    gaba_psc->reset();
}
void Neuron::print_g(ostream &s)
{
    ampa_psc->print_state(s);
    nmda_psc->print_state(s);
    gaba_psc->print_state(s);

    for(int i=0; i<currentstep; i++)
        s<<gi[i+1]<<'\t';
    s<<endl;

    for(int i=0; i<currentstep; i++)
        s<<ge[i+1]<<'\t';
    s<<endl;
}
void Neuron::print_voltage(ostream &s)
{

    for(int i=0; i<currentstep; i++)
        s<<voltage[i+1]<<'\t';
    s<<endl;
}
void Neuron::print_spiketime(ostream &s)
{
    if(spike_time.size()<=0)
        return;
    for(int i=0; i<spike_time.size(); i++)
        s<<spike_time[i]<<endl;
    s<<endl;
}

#endif
