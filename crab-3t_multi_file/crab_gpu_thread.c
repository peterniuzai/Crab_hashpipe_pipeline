/*crab_gpu_thread.c
 *
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include "hashpipe.h"
#include "crab_databuf.h"
#include "math.h"
static void *run(hashpipe_thread_args_t * args)
{
    // Local aliases to shorten access to args fields
    crab_input_databuf_t *db_in = (crab_input_databuf_t *)args->ibuf;
    crab_output_databuf_t *db_out = (crab_output_databuf_t *)args->obuf;
    hashpipe_status_t st = args->st;
    const char * status_key = args->thread_desc->skey;

    int rv;
    uint64_t mcnt=0;
    int curblock_in=0;
    int curblock_out=0;

//	short int xx_tmp[N_CHANS_SPEC];
//	short int yy_tmp[N_CHANS_SPEC];
//	short int xx,yy;
	short int I[N_CHANS_BUFF/N_POLS_CHAN];
	//short int xy_real_tmp[N_CHANS_SPEC];
	//short int xy_img_tmp[N_CHANS_SPEC];
//	short int packet_tmp[N_CHANS_SPEC*N_SPEC_BUFF*N_POLS_CHAN];

/*
	I = X*X + Y*Y
	Q = X*X - Y*Y
	U = 2 * XY_real
	V = -2 * XY_img
*/
/*	struct Full_Stokes{
		short int I[N_CHANS_SPEC];
		short int Q[N_CHANS_SPEC];
		short int U[N_CHANS_SPEC];
		short int V[N_CHANS_SPEC]; 
	}*/
    while (run_threads()) {

        hashpipe_status_lock_safe(&st);
        hputi4(st.buf, "COVT-IN", curblock_in);
        hputs(st.buf, status_key, "waiting");
        hputi4(st.buf, "COVT-OUT", curblock_out);
	hputi8(st.buf,"COV-MCNT",db_in->block[curblock_in].header.mcnt);
        hashpipe_status_unlock_safe(&st);

        // Wait for new input block to be filled
        while ((rv=crab_input_databuf_wait_filled(db_in, curblock_in)) != HASHPIPE_OK) {
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

        // Got a new data block, update status and determine how to handle it
        /*hashpipe_status_lock_safe(&st);
        hputu8(st.buf, "GPUMCNT", db_in->block[curblock_in].header.mcnt);
        hashpipe_status_unlock_safe(&st);*/


        // Note processing status
        hashpipe_status_lock_safe(&st);
        hputs(st.buf, status_key, "processing");
        hashpipe_status_unlock_safe(&st);
	for(int j=0;j<N_CHANS_BUFF/N_POLS_CHAN;j++){
		I[j]=sqrt(pow(db_in->block[curblock_in].data[j*N_POLS_CHAN],2) + pow(db_in->block[curblock_in].data[j*N_POLS_CHAN+1],2));
                }

        // Mark input block as free and advance
        crab_input_databuf_set_free(db_in, curblock_in);
        curblock_in = (curblock_in + 1) % db_in->header.n_block;

	// Wait for new output block to be free
	 while ((rv=crab_output_databuf_wait_free(db_out, curblock_out)) != HASHPIPE_OK) {
            if (rv==HASHPIPE_TIMEOUT) {
                hashpipe_status_lock_safe(&st);
                hputs(st.buf, status_key, "blocked out");
                hashpipe_status_unlock_safe(&st);
                continue;
            } else {
                hashpipe_error(__FUNCTION__, "error waiting for free databuf");
                pthread_exit(NULL);
                break;
            }
        }
	memcpy(db_out->block[curblock_out].I,I,N_CHANS_BUFF/N_POLS_CHAN*sizeof(short int));
        // Mark output block as full and advance
        crab_output_databuf_set_filled(db_out, curblock_out);
        curblock_out = (curblock_out + 1) % db_out->header.n_block;
		
	mcnt++;
        /* Check for cancel */
        pthread_testcancel();
    }
    return THREAD_OK;
}

static hashpipe_thread_desc_t crab_gpu_thread = {
    name: "crab_gpu_thread",
    skey: "COV-STAT",
    init: NULL,
    run:  run,
    ibuf_desc: {crab_input_databuf_create},
    obuf_desc: {crab_output_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&crab_gpu_thread);
}

