/*
 * crab_output_thread.c
 * 
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include "hashpipe.h"
#include "crab_databuf.h"
#include "filterbank.h"
#include <sys/time.h>
//#include "crab_net_thread.h"
//#include "crab_net_thread.c"
extern unsigned long long  miss_gap;
//printf
//extern bool store_flag
static void *run(hashpipe_thread_args_t * args)
{
	//printf("\n%d Mbytes for each Filterbank file.\n ",N_MBYTES_PER_FILE);
	// Local aliases to shorten access to args fields
	extern unsigned long long  miss_pkt1;
	// Our input buffer happens to be a crab_ouput_databuf
	crab_output_databuf_t *db = (crab_output_databuf_t *)args->ibuf;
	hashpipe_status_t st = args->st;
	const char * status_key = args->thread_desc->skey;
	int rv;
	int block_idx = 0;
	long long unsigned N_Mbytes_save = 0;
	long long unsigned N_Mbytes_file = N_MBYTES_PER_FILE;
	//double N_Mbytes_save = 0;
	//double N_Mbytes_file = N_MBYTES_PER_FILE;
	int filb_flag = 1;
	FILE * crab_file;
/*
	I = X*X + Y*Y
	Q = X*X - Y*Y
	U = 2 * XY_real
	V = -2 * XY_img
*/
//	short int I[N_CHANS_SPEC*N_SPEC_BUFF];
	//short int Q[N_CHANS_SPEC];
	//short int U[N_CHANS_SPEC];
	//short int V[N_CHANS_SPEC]; 
/*	char	f_fil[256];
	struct tm  *now;
	time_t rawtime;
	FILE * crab_file;
	printf("\n\nopen new filterbank file...\n\n");
	time(&rawtime);
	now = localtime(&rawtime);
	strftime(f_fil,sizeof(f_fil), "./data_%Y-%m-%d_%H-%M-%S.fil",now);
	WriteHeader(f_fil);
	printf("write header done!\n");
	crab_file=fopen(f_fil,"a+");
	printf("starting write data...\n");	
*/
	/* Main loop */
	while (run_threads()) {

		hashpipe_status_lock_safe(&st);
		hputi4(st.buf, "OUTBLKIN", block_idx);
		hputi8(st.buf, "DATSAVMB",N_Mbytes_save);
		hputs(st.buf, status_key, "waiting");
		hashpipe_status_unlock_safe(&st);

		// get new data
		while ((rv=crab_output_databuf_wait_filled(db, block_idx))
		!= HASHPIPE_OK) {
		if (rv==HASHPIPE_TIMEOUT) {
			hashpipe_status_lock_safe(&st);
			hputs(st.buf, status_key, "blocked");
			hashpipe_status_unlock_safe(&st);
			continue;
			} else {
				hashpipe_error(__FUNCTION__, "error waiting for filled databuf");
				pthread_exit(NULL);
				break;
			}
		}

		hashpipe_status_lock_safe(&st);
		hputs(st.buf, status_key, "processing");
		hashpipe_status_unlock_safe(&st);

		//memcpy(Full_Stokes,db->block[block_idx].Stokes_Full,PKTSIZE*sizeof(char));
//		memcpy(I,db->block[block_idx].I,N_CHANS_SPEC*N_SPEC_BUFF*sizeof(short int));
		//memcpy(Q,db->block[block_idx].Q,N_CHANS_SPEC*sizeof(short int));
		//memcpy(U,db->block[block_idx].U_REAL,N_CHANS_SPEC*sizeof(short int));
		//memcpy(V,db->block[block_idx].V,N_CHANS_SPEC*sizeof(short int));
		//printf("first data from output thread is:%c\n",I[0]);

		if (miss_gap >= 2){
		//	printf("\nready to write\n");
			if (filb_flag ==1){
			        char    f_fil[256];
			        struct tm  *now;
				time_t rawtime;

			        printf("\n\nopen new filterbank file...\n\n");
	        		time(&rawtime);
			        now = localtime(&rawtime);
		        	strftime(f_fil,sizeof(f_fil), "./data_%Y-%m-%d_%H-%M-%S.fil",now);
	        		WriteHeader(f_fil);
			        printf("write header done!\n");
			        crab_file=fopen(f_fil,"a+");
		        	printf("starting write data to %s...\n",f_fil);
						}
	
			fwrite(db->block[block_idx].I,sizeof(short int),N_CHANS_BUFF/N_POLS_CHAN,crab_file);
			N_Mbytes_save += BUFF_SIZE/N_POLS_CHAN/1024/1024;
			
	//		printf("\ndata save:%lld\n",N_Mbytes_save);
	//		printf("\nfile save:%lld\n",N_Mbytes_file);
	//		printf("\ndevide?:%lld\n",N_Mbytes_save%N_Mbytes_file);
			//if (N_Mbytes_save % N_Mbytes_file ==0){
			//		filb_flag = 1;
			//		}
			//else{
			//		filb_flag = 0;
			//		}		
			filb_flag = 0;
		
				}	
		//fwrite(Q,sizeof(short int),N_CHANS_SPEC,crab_file);
		//fwrite(U,sizeof(short int),N_CHANS_SPEC,crab_file);
		//fwrite(V,sizeof(short int),N_CHANS_SPEC,crab_file);
		//char a = Full_Stokes[0];
		//fwrite(&a,sizeof(a),1,crab_file);
		//fwrite(&Full_Stokes[0],sizeof(a),1,crab_file);
		//sleep(0.1);
		//hashpipe_status_lock_safe(&st);
		//hputi4(st.buf, "OUT_FIRST", I[0]);
		//hashpipe_status_unlock_safe(&st);

		crab_output_databuf_set_free(db,block_idx);
		block_idx = (block_idx + 1) % db->header.n_block;
		

		//Will exit if thread has been cancelled
		pthread_testcancel();

	}
	fclose(crab_file);
	return THREAD_OK;
}

static hashpipe_thread_desc_t crab_output_thread = {
	name: "crab_output_thread",
	skey: "OUTSTAT",
	init: NULL, 
	run:  run,
	ibuf_desc: {crab_output_databuf_create},
	obuf_desc: {NULL}
};

static __attribute__((constructor)) void ctor()
{
	register_hashpipe_thread(&crab_output_thread);
}

