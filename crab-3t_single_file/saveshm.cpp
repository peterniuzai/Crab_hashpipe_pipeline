//g++ filterbank.cpp saveshm.cpp -o filterbank
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "filterbank.h"
#include <sys/time.h>

int main(void)
{
	double Fch1=1699.609375;
	double Foff=-0.78125;
	int Nchans=384;
	float Tsamp=64e-6;
	int Nbits=8;
	int Nifs=4;
	int start=time((time_t*)NULL);

	/*double Year, Month, Day;
	double Mjd;
	time_t timep;
	struct tm *p;
	time(&timep);
	p=gmtime(&timep);
	Year=p->tm_year+1900;
	Month=p->tm_mon+1;
	Day=p->tm_mday;
	Mjd = UTC2MJD(Year, Month, Day); 
	gettimeofday( &currenttime, NULL );
	time(&timep);
	p=gmtime(&timep);
	MJD=Mjd+(double)((p->tm_hour)/24.0)
                               +(double)(p->tm_min/1440.0)
                               +(double)(p->tm_sec/86400.0)
                               +(double)(currenttime.tv_usec/86400.0/1000000.0);

	printf("MJD time is:%.15f\n",(ph+i)->MJD);*/

	const char * fname= "filterbank-test";
	//Start to allocate the header information
	//long int header_size=sizeof(_MEMELEM)*Nseg;
	FilterBankData fil;
	strcpy(fil.Source_name, "J0000-0000");
	fil.UseFrequencyTable=false;
	//fil.Telescope_id=1;
	//fil.Machine_id=1;
	fil.Data_type=1;
	fil.Az_start=0;
	fil.Za_start=0;
	fil.Src_dej=0;
	fil.Src_raj=0;
	fil.Tstart=50011.12445566;
	fil.Tsamp=Tsamp;
	fil.Fch1=Fch1;
	fil.Foff=Foff;
	fil.Nchans=Nchans;
	fil.Nbits=Nbits;
	fil.Nbeams=1;
	fil.Nifs=Nifs;
	fil.RefDM=0;
	//fil.Nsamples=Nsl/Nchans/Nifs;
	//fil.pData=(float *)pd_proc;
	char outputname[4096];
	sprintf(outputname, "%s.fil", fname);
	fil.WriteToFile(outputname);
	
	usleep(100);
	exit(0);
}


