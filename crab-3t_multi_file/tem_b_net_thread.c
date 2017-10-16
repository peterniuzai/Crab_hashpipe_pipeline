/*
 * crab_net_thread.c
 *
 *  
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include "hashpipe.h"
#include "crab_databuf.h"

//defining a struct of type hashpipe_udp_params as defined in hashpipe_udp.h
static struct hashpipe_udp_params params;

static int init(hashpipe_thread_args_t * args)
{
        hashpipe_status_t st = args->st;
        //strcpy(params.bindhost,"127.0.1.1");
	params.bindport = 10000;
	strcpy(params.bindhost,"10.10.12.35");
        //selecting a port to listen to
        //params.bindport = 10000;
        //strcpy(params.bindhost,"10.10.12.35");
        //selecting a port to listen to
        //params.bindport = 10000;
        //params.packet_size = 0;
        hashpipe_udp_init(&params);
        hashpipe_status_lock_safe(&st);
        hputi8(st.buf, "NPACKETS", 0);
        hputi8(st.buf, "NBYTES", 0);
        hashpipe_status_unlock_safe(&st);
        return 0;

}


static void *run(hashpipe_thread_args_t * args)
{
    crab_input_databuf_t *db  = (crab_input_databuf_t *)args->obuf;
    hashpipe_status_t st = args->st;
    const char * status_key = args->thread_desc->skey;

    /* Main loop */
    int i, rv,input,n;
    //uint64_t mcnt = 0;
    unsigned long long mcnt	= 0;
    unsigned long long seq	= 1;
    unsigned long long mcnt_i    = 0;
    double rcvmb    = 0;
    int	start_flag  = 0;
    double miss_rate = 0;
    int block_idx   = 0;
    unsigned long long miss_pkt    = 0;
    int ic =0;//temp variable
    char *data0;
    data0 = (char *)malloc(PKTSIZE*sizeof(char));
    unsigned long long  npackets = 0; //number of received packets
    unsigned long long  nbytes = 0;  //number of received bytes
    

    while (run_threads()) {

	//rcvmb  = nbytes/1024/1024;
        hashpipe_status_lock_safe(&st);
        hputs(st.buf, status_key, "waiting");
        hputi4(st.buf, "NETBKOUT", block_idx);
//	hputi8(st.buf,"NETMCNT",mcnt);
        hputi8(st.buf, "NPACKETS", npackets);
        hputi8(st.buf, "RCVMB", rcvmb);
        hashpipe_status_unlock_safe(&st);
 
        // Wait for data
        /* Wait for new block to be free, then clear it
         * if necessary and fill its header with new values.
         */
        while ((rv=crab_input_databuf_wait_free(db, block_idx)) 
                != HASHPIPE_OK) {
            if (rv==HASHPIPE_TIMEOUT) {
                hashpipe_status_lock_safe(&st);
                hputs(st.buf, status_key, "blocked");
                hashpipe_status_unlock_safe(&st);
                continue;
            } else {
                hashpipe_error(__FUNCTION__, "error waiting for free databuf");
                pthread_exit(NULL);
                break;
            }
        }

        hashpipe_status_lock_safe(&st);
        hputs(st.buf, status_key, "receiving");
        hashpipe_status_unlock_safe(&st);

        n = recvfrom(params.sock,data0,PKTSIZE*sizeof(char),0,NULL,NULL);
	//fprintf(stderr,"\nStarting to receiving data...\n\n");
	//if received packet has data,count packtes,mark it as the first or second packet. 
	if (n>0){
		
//        	printf("received %d bytes from number %ld packets!\n",n,mcnt);
        	//sleep(0.1);
		memcpy(&mcnt,data0,N_BYTES_COUNTER*sizeof(char));
		//fprintf(stderr,"1st print:mcnt number :%lli\n, seq number:%lli\n\n",mcnt,seq);
		seq =  (mcnt /512)%8;
		
		//if (ic == 16){
		//	exit(1);
		//	}else{
		//	ic++;	}

		fprintf(stderr,"mcnt number :%lli\n, seq number:%lli\n,miss_p:%lld,total:%lld ,mis_ra:%5.4f\n",mcnt,seq,miss_pkt,npackets,miss_rate);
		//fprintf(stderr,"first seq number :%lli\n\n",seq);
		if (seq == 0 || start_flag == 1){
			nbytes += n;
			if (npackets=0){mcnt_i = mcnt;}

			if ((mcnt_i)  != mcnt){
				loss_p = (mcnt - mcnt_i)/N_CHAN_PER_PACK;
				loc    = mcnt_i/512%8 ;
				if (loc+loss_p <8){
					memcpy(db->block[block_idx].data+(loc*DATA_SIZE_PACK),0,loss_p*DATA_SIZE_PACK*sizeof(char));//problem
						}
				else{
					
					memcpy(db->block[block_idx].data+(loc*DATA_SIZE_PACK),0,loss_p*DATA_SIZE_PACK*sizeof(char));}
				mcnt_i = mcnt;
                        	miss_pkt +=1;
				miss_rate = (double)miss_pkt/npackets;
					 }
			else{
				//fprintf(stderr,"received %d bytes, sequency number: %lli!\n",n,seq);
				//exit(1);
				memcpy(db->block[block_idx].data+seq*DATA_SIZE_PACK,data0+8,DATA_SIZE_PACK*sizeof(char));
				}
			db->block[block_idx].header.mcnt = mcnt_i;
	        	hashpipe_status_lock_safe(&st);
			hputi8(st.buf,"NETMCNT",mcnt);
			hputi8(st.buf,"PKTseq",seq);
			hputi8(st.buf,"MiSSPKT",miss_pkt);
			hputi8(st.buf,"MiSSPKT%",miss_rate*1000000);
		        hashpipe_status_unlock_safe(&st);

        		// Mark block as full
	        	if(crab_input_databuf_set_filled(db, block_idx) != HASHPIPE_OK) {
                	hashpipe_error(__FUNCTION__, "error waiting for databuf filled call");
	                pthread_exit(NULL);
						    				        }
			start_flag = 1;
			mcnt_i +=512;
			npackets++;
						}
		}else{continue;}

        // Setup for next block
        block_idx = (block_idx + 1) % db->header.n_block;
//        mcnt++;

        /* Will exit if thread has been cancelled */
        pthread_testcancel();
    }

    // Thread success!
    return THREAD_OK;
}

static hashpipe_thread_desc_t crab_net_thread = {
    name: "crab_net_thread",
    skey: "NETSTAT",
    init: init,
    run:  run,
    ibuf_desc: {NULL},
    obuf_desc: {crab_input_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&crab_net_thread);
}
