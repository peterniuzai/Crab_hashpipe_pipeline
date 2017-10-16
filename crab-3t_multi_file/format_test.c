#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <byteswap.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "filterbank.h"
#include <sys/time.h>
#define _CHAR_SWAP_SIZE 256

int main(void)
{
	
}
#define PKTSIZE 4104
void main(void){
FILE * crab_file,fil_file;
short int a;
char raw_paket[PKTSIZE];
short int raw_data[2048];
long int counter;
short int xx[512];
short int yy[512];
short int xy_real[512];
short int xy_img[512];
crab_file=fopen("./data_2017-06-06_16-36-36","r");
WriteHeader("format_test");

fil_file=fopen("./format_test.fil","a+");
//fwrite(&a,sizeof(a),1,crab_file);
for(int i=0;i<10;i++){
	//fread(&counter,sizeof(long int),1,crab_file);
	fread(raw_data,sizeof(short int),2048,crab_file);
	for(int j=0;j<512;j++){
		xx[512-j-1]=raw_data[j*4];
		yy[512-j-1]=raw_data[j*4+1];
		xy_real[512-j-1]=raw_data[j*4+2];
		xy_img[512-j-1]=raw_data[j*4+3];
	}
	//counter = __bswap_64(counter);
	//counter = counter & 0xffff;
	fwrite(xx,sizeof(short int),512,fil_file);
	fwrite(yy,sizeof(short int),512,fil_file);
	fwrite(xy_real,sizeof(short int),512,fil_file);
	fwrite(xy_img,sizeof(short int),512,fil_file);
	printf("raw_data is:%d\n",raw_data[0]);
	printf("size of counter is: %ld\n",sizeof(raw_data[0]));
	}
}
