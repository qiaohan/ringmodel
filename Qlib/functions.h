#ifndef FUN_H
#define FUN_H

double regular_s(double s0,double h,double tao)
{ double s1;
  s1=s0*exp(-h/tao);
  return s1;
}

/*
double regular_g(double s0,double g0,double h,double tao)
{double g1=g0*exp(-h/tao)+s0/tao*h*exp(-h/tao);
 return g1;
}
*/
double regular_g2(double s0,double g0,double h,double tao1,double tao2)
{
  double g1=g0*exp(-h/tao2)+s0*tao1/(tao1-tao2)*(exp(-h/tao1)-exp(-h/tao2));
  return g1;
 }

double add_s(double s0,double h,double ts,double str,double tao)
{ double s1;
  s1=s0+str/tao*exp(-(h-ts)/tao);
  return s1;
}
/*
double add_g(double g0,double h,double ts,double str,double tao)
{double g1;
  g1=g0+str/(tao*tao)*(h-ts)*exp(-(h-ts)/tao);
  return g1;
}
*/
double add_g2(double g0,double h,double ts,double str,double tao1,double tao2)
{
	double g1=g0+str/(tao1-tao2)*(exp(-(h-ts)/tao1)-exp(-(h-ts)/tao2));
	return g1;
}

double find_ts(double vthr,double v1,double v2,double h)
{double ts;
 ts=(vthr-v1)/(v2-v1)*h;
 return ts;
}

double find_v2(double g_e1,double g_e2,double g_i1,double g_i2,double ts,double h)
{double v1,v2,k1,k2;
 double alfa1,alfa2,beta1,beta2;
 alfa1=50+g_e1+g_i1;
 alfa2=50+g_e2+g_i2;
 beta1=(double)14/3*g_e1-(double)2/3*g_i1;
 beta2=(double)14/3*g_e2-(double)2/3*g_i2;

 v1=-ts*(beta1+beta2-alfa2*beta1*h)/2/
    (1+ts*(-alfa1-alfa2+alfa1*alfa2*h)/2);
 k1=-alfa1*v1+beta1;
 k2=-alfa2*(v1+k1*h)+beta2;
 v2=v1+h*(k1+k2)/2;
 //if(v2>1)
 //	 v2=0;
 return v2;
}

double RK2(double v1,double ge1,double ge2,double gi1,double gi2,double h)
{double v2;
 double k1,k2;
 k1=-(50+ge1+gi1)*v1+(double)14/3*ge1-(double)2/3*gi1;
 k2=-(50+ge2+gi2)*(v1+h*k1)+(double)14/3*ge2-(double)2/3*gi2;
 v2=v1+(k1+k2)*h/2;

 return v2;
}
/*
void write_data_to_file(const double* array, int arrayLen, const char* fileName)
{
FILE* fp;
int i;
  char outputfilepath[100];
  strcpy(outputfilepath,outputpath);
fp = fopen(strcat(outputfilepath,fileName), "w");
if (fp == NULL)
return;

for (i=0; i<arrayLen; i++)
fprintf(fp, "%f\r\n", array[i]);
fclose(fp);
}

void write_data_to_file1(const double** array, int row, int column, const char* fileName,int blank)
{
  FILE* fp;
  int i, j;
   char outputfilepath[100];
  strcpy(outputfilepath,outputpath);
  fp = fopen(strcat(outputfilepath,fileName), "w");
  if (fp == NULL)
    return;
 
  for (i=0; i<row; i++)
  {for (j=0; j<column-1; j++)
   {//if(j%blank==0&&j!=0)
	if(j%blank==0)
    {fprintf(fp, "%10.6f ", *((double*)array+i*column+j));}
    }
    fprintf(fp, "\n");
   }
 
  fclose(fp);
}

void write_data_to_file2(const double** array, int row, int column, const char* fileName)
{
  FILE* fp;
  int i, j;
  char outputfilepath[100];
  strcpy(outputfilepath,outputpath);
  fp = fopen(strcat(outputfilepath,fileName), "w");
  if (fp == NULL)
    return;
 
  for (i=0; i<row; i++)
  {for (j=0; j<column; j++)
   {fprintf(fp, "%10.9f ", *((double*)array+i*column+j));}
    fprintf(fp, "\n");
   }
 
  fclose(fp);
}
*/
char* itoa(int val, int base=10){
	
	static char buf[32] = {0};
	
	int i = 30;
	
	for(; val && i ; --i, val /= base)
	
		buf[i] = "0123456789abcdef"[val % base];
	
	return &buf[i+1];
	
}

#endif
