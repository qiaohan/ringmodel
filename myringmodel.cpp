#include "Qlib/neuron_class.h"
#include "Qlib/connection_class.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <stdlib.h>
using namespace std;


int main(int argc, char* argv[])
{
    /*
     * Init neurons independently
     */
    ofstream f;
    f.open("neuron_parameter.txt");
    Neuron * neurons[NUM_NEURON];
    srand( (int)time(0) );
    for(int i=0; i<NUM_NEURON; i++)
    {
        //int excit = rand()%1000<750? 1:0;
        int excit = i%4!=0?1:0;
        int sim = rand()%1000<500? 1:0;
        double phase = (double)( rand()%1000 )/1000.0 *2*PI;
        double theta = PI*((double)i)/NUM_NEURON;
        neurons[i] = new Neuron(phase,theta,sim,excit);
        neurons[i]->initvoltage(0.0);
        neurons[i]->initpsc("ampa",0,0);
        neurons[i]->initpsc("nmda",0,0);
        neurons[i]->initpsc("gaba",0,0);
        f<<phase<<'\t'<<theta<<'\t'<<sim<<'\t'<<excit<<endl;
    }
    f.close();
    /*
     * init connection
     */
    
    Connection con(NUM_NEURON,neurons,SPARSE);
    con.saveConnect("connect.txt");
    /*
     *loop for evoluation:
      1.update g w.r.t. lgn spike
      2.update v w.r.t g, and find the spike
      3.update g w.r.t v's spike
      But the v is not updated in this loop, it will cause residual for v and v's spike
      Luckily, accurate g is enough for next loop, and residual for v overall during the whole evolution
      will be order cubic time step.
     * */
    vector<double> neuspiketime[NUM_NEURON];
    for(int cycle=0; cycle<NUM_CYCLE; cycle++)
    {
         for(int sec=0; sec<SECOND_PER_CYCLE; sec++)
         {
            int second = sec + cycle*SECOND_PER_CYCLE;
            cout<<"time: "<<second<<" s"<<endl;
            for(int t=0; t<CLK; t++)
            {
                 double time = second + t*TINVAL;
                 double simulate_theta = PI*sec/SECOND_PER_CYCLE;
                 double simulate_phase = w*time;
                 double simulate_constrast = 12;
                 for(int i=0; i<NUM_NEURON; i++)
                 {
                     neurons[i]->evolve(simulate_theta,simulate_phase,simulate_constrast);
                 }
                /*
                *   check the spike for every neuron(reset the voltage and refraction) 
                *   and convey the spike to connected others
                * */
                for(int i=0; i<NUM_NEURON; i++)
                 {
                     double tts;
                        if(neurons[i]->checkspike(&tts))
                        {
                            //cout<<time<<'\t';
                            //convey to others
                            neuspiketime[i].push_back(time+tts-TINVAL); 
                            for(int j = 0; j<NUM_NEURON; j++)
                            {
                                if( neurons[i]->isexcit() )
                                    neurons[j]->Espikecome(0*con.Estrong[i][j],tts);
                                else
                                    neurons[j]->Ispikecome(0*con.Istrong[i][j],tts);
                            }
                        }
                 }
            }

             
             
             /*
             * one second finished, store the states, cal and store the firing rates 
             * and reset the neurons for next second evolving
             * */

            for(int i =0; i<NUM_NEURON; i++)
            {
                 string s1 = itoa(i+1);
                 string s2 = itoa(second+1);
                 string fname(argv[1]);
                 fname = fname + "/voltage_neu_"+s1+"sec_"+s2+".txt";
                 /*
                 ofstream f;
                 f.open(fname.c_str());
                 neurons[i]->print_voltage(f);
                 f.close();
                 */
                 neurons[i]->reset();
            }
         }
    }

    cout<<"evolve finished"<<endl;
            for(int i =0; i<NUM_NEURON; i++)
            {
                 string s1 = itoa(i+1);
                 string fname(argv[1]);
                 fname = fname +  "/spike_time_neu_"+s1+".txt";
                 ofstream f;
                 f.open(fname.c_str());
                 //neurons[i]->print_spiketime(f);
                 for(int k = 0; k<neuspiketime[i].size(); k++)
                     f<<neuspiketime[i][k]<<endl;
                 f.close();
            }
    
    cout<<"all finished!"<<endl;
    for(int i=0; i<NUM_NEURON; i++)
    {
        delete neurons[i];
    }

    return 0;
}
