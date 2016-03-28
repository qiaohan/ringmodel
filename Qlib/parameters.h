#ifndef PARAMETER_H
#define PARAMETER_H

#define pi 3.1415926535897932384626433832795
#define PI 3.1415926535897932384626433832795

double w=50.2655;//时间频率=8

#define AMPAR 0.001
#define AMPAD 0.005
#define NMDAR 0.002
#define NMDAD 0.08
#define GABAR 0.001
#define GABAD 0.01

#define ALPHA 1

#define REF_E 0.003
#define REF_I 0.001

#define VTHR 1
//voltage threshold

#define NUM_CYCLE 1
#define SECOND_PER_CYCLE 8
#define CLK 10000 
#define TINVAL 0.0001 
//step time
#define NUM_NEURON 1024

 //double sto_str=0.28;//LGN连接强度
 double lgnE_str=0.1;
 double lgnI_str=1;//随机抑制电导连接强度
 //double f_sinh=190;//随机抑制电导放电率
 double lgnIfr=180;

 //double S_EE0=2;double S_EE1=6.8;//简单，复杂细胞的兴奋性连接强度
 //double S_EI0=4.5;double S_EI1=9.5;
 double S_EE0=2;double S_EE1=7.4;
 double S_EI0=4.5;double S_EI1=9.5;
 //double S_E[1024]={0};//整体细胞的兴奋性连接强度
 //double S_I=8;
 //double S_I=8;
#define SI 8
 //double sparse=0.06;//定义偶联疏密度
 //double sparse=0.08;
 //double Success=1;//突触链接的成功SP率
#define SPARSE 0.08
#define SUCCESS 1

#endif
