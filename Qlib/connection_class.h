#ifndef CONNECTION_H
#define CONNECTION_H

#include "neuron_class.h"
#include <fstream>

class Connection
{
    public:
        Connection(int num, Neuron ** neu, double sparse);
        ~Connection();
        // 2-d strong array, Istrong[i][j] is the inhibitory degree from neuron i to j
        double ** Istrong;
        double ** Estrong;
        void saveConnect(string fname);
    private:
        int number;
};
void Connection::saveConnect(string fname)
{
    ofstream f;
    f.open(fname.c_str());
    for(int i=0;i<number;i++)
    {
        if(Estrong[i])
        {
            for(int j = 0; j<number; j++)
                f<<Estrong[i][j]<<'\t';
        }
        else
        {
            for(int j=0; j<number; j++)
                f<<(-Istrong[i][j])<<'\t';
        }
        f<<endl;
    }
    f.close();
}
Connection::Connection(int num, Neuron ** neu, double sparse)
{
    number = num;
    Estrong = new double*[num];
    Istrong = new double*[num];
    double S_E,S_I=SI;
    for(int i=0; i<num; i++)
    {
        if( neu[i]->isexcit() )
        {
            Estrong[i] = new double[num];
            Istrong[i] = NULL;
            for(int j=0; j<num; j++)
            { 
                if(neu[j]->isexcit())
                    S_E = (1-neu[j]->issimple())*S_EE1+S_EE0;
                else
                    S_E = (1-neu[j]->issimple())*S_EI1+S_EI0;
                //excitatory Guassion, 3*num/4 neurons
                double sigma = 220;
                int d = j-i;
                d = d>0?d:(-d);
                int d1 = d;
                int d2 = NUM_NEURON - d;
                int neu_distance = d1<d2?d1:d2;
                double musquare = neu_distance*neu_distance;
                double ke = 1/(sigma*sqrt(2*PI))*exp(-musquare/(2*sigma*sigma))*4.0/3.0;
            
                Estrong[i][j] = ke*S_E;
                int p = rand()%1000;
                if(p<sparse*1000)
                    p=1;
                else
                    p=0;
            
                //Estrong[i][j] *= p;
                //Estrong[i][j] /= sparse;
            }
        }
        else
        {
            Istrong[i] = new double[num];
            Estrong[i] = NULL;
            for(int j=0; j<num; j++)
            { 
                double ki = 4.0/num; 
                Istrong[i][j] = ki*S_I;
                int p = rand()%1000;
                if(p<sparse*1000)
                    p=1;
                else
                    p=0;
            
                //Istrong[i][j] *= p;
                //Istrong[i][j] /= sparse;
            }
        }
    }
}
Connection::~Connection()
{
    for(int i=0; i<number; i++)
    {
        if(Estrong[i])
            delete [] Estrong[i];
        if(Istrong[i])
            delete [] Istrong[i];
    }
    delete [] Estrong;
    delete [] Istrong;
}

#endif

