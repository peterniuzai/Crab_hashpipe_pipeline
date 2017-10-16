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
//#include "crab_net_thread.h"
//defining a struct of type hashpipe_udp_params as defined in hashpipe_udp.h
static struct hashpipe_udp_params params;
//unsigned long long miss_pkt = 0;
long miss_gap = 0;
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
	hputi8(st.buf,"NETMCNT",0);
        hputi8(st.buf, "NPACKETS", 0);
        hputi8(st.buf, "RCVMB",0);
	hputi8(st.buf, "COV-MCNT",0);
	hputi8(st.buf,"MiSSPKT",0);
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
    unsigned long long mcnt	 = 0;
    unsigned long long mcnt_i    = 0;
    unsigned long long rcvmb  	 = 0;
    unsigned long long nbytes	 = 0;  //number of received bytes
    long npackets		 = 0;  //number of received packets
    long miss_pkt	= 0;
    long offset		= 0;
 //   double miss_rate	= 0;
    int start_flag	= 0;
    int miss_spec	= 0;
    int block_idx	= 0;
    int seq		= 0;
    bool store_flag	= 0;
    //int ic =0;//temp variable
    char *data0;
    data0 = (char *)malloc(PKTSIZE*sizeof(char));
    //sleep(5);

    while (run_threads()) {

	rcvmb  = nbytes/1024/1024;
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
	if (n == PKTSIZE){
        	//sleep(0.1);
		memcpy(&mcnt,data0,N_BYTES_COUNTER*sizeof(char));
		seq =  (mcnt /N_CHAN_PER_PACK)%N_PACKETS_PER_SPEC;
		//printf("\n###start!##\nseq number:%lli,miss_p:%lld,npackets:%lld ,mcnt:%lld,mcnt_i:%lld,flag:%d\n",seq,miss_pkt,npackets,mcnt,mcnt_i,start_flag);
		if (seq == 0 || start_flag == 1){
			if (npackets == 0 ){mcnt_i = mcnt;}
			if ( miss_gap == 1 && seq == 0 ){mcnt_i = mcnt;offset=0;miss_gap++;fprintf(stderr,"\n\nhahaha!!!! miss_gap:%ld\n\n",miss_gap);}//napackets always start with 0 it works. problem?
			nbytes += n;
			npackets++;
			if (mcnt_i < mcnt){
				if (miss_gap >=2){
					fprintf(stderr,"heloo??????,mcnt%lld,mcnt_i:%lld,miss_pkt:%ld,seq:%d \n",mcnt,mcnt_i,miss_pkt,seq);
					//exit(1);
					miss_pkt +=(mcnt-mcnt_i)/N_CHAN_PER_PACK;
					fprintf(stderr,"misspkt:%ld\noffset:%ld\n",miss_pkt,offset);
					//exit(1);
					if ((offset + miss_pkt * N_CHAN_PER_PACK * N_POLS_CHAN) >= N_CHANS_BUFF){
						printf("\nProblem?\n");
						//exit(1);
						memset(db->block[block_idx].data+offset,0,(N_CHANS_BUFF-offset)*N_BYTES_DATA_POINT*sizeof(char));
                                                db->block[block_idx].header.mcnt = mcnt_i;
                                                offset  = 0;    
						//start_flag = 0;
                                                mcnt_i += (N_CHANS_BUFF-offset)/N_POLS_CHAN;
						// Mark block as full
                                                if(crab_input_databuf_set_filled(db, block_idx) != HASHPIPE_OK) {
                                                        hashpipe_error(__FUNCTION__, "error waiting for databuf filled call");
                                                        pthread_exit(NULL);}
                                                        block_idx = (block_idx + 1) % db->header.n_block;
					}else{
 						 memset(db->block[block_idx].data+offset,0, miss_pkt*DATA_SIZE_PACK *sizeof(char));
						 printf("Oh no!");
						 mcnt_i = mcnt+N_CHAN_PER_PACK;
						 //exit(1);
						 db->block[block_idx].header.mcnt = mcnt_i;
						 offset += miss_pkt * N_CHAN_PER_PACK * N_POLS_CHAN;
						}
							}
				miss_gap++;
				printf("####How is miss_gap:%ld#####\n",miss_gap);
				if (miss_gap == 1){start_flag = 0;}
				//exit(1);
				//miss_rate = (double)miss_pkt/npackets;
						}//if(mcnt_i !=mcnt)
			else{  
				//fprintf(stderr,"Processing no loos packet..\n\nmcnt:%lld,mcnt_i:%lld,npackets:%lld,seq:%lld\n",mcnt,mcnt_i,npackets,seq);
				
				memcpy(db->block[block_idx].data+offset,data0+N_BYTES_COUNTER,DATA_SIZE_PACK*sizeof(char));
				db->block[block_idx].header.mcnt = mcnt_i;
				hashpipe_status_lock_safe(&st);
				hputi8(st.buf,"NETMCNT",mcnt);
				hputi8(st.buf,"PKTseq",seq);
				hputi8(st.buf,"MiSSPKT",miss_pkt);
				//hputi8(st.buf,"MiSSPKT%",miss_rate*1000000);
				hashpipe_status_unlock_safe(&st);
			
				offset += DATA_SIZE_PACK/N_BYTES_DATA_POINT;
				start_flag = 1;
        	                mcnt_i += N_CHAN_PER_PACK;
				
				


				if (offset == N_CHANS_BUFF ){
                                // Mark block as full
                                        if(crab_input_databuf_set_filled(db, block_idx) != HASHPIPE_OK) {
                                                hashpipe_error(__FUNCTION__, "error waiting for databuf filled call");
                                                pthread_exit(NULL);
                                                                                                       }
                                        block_idx = (block_idx + 1) % db->header.n_block;
					//printf("\n\nWerid!!!\n\n");
                                        //exit(1);
                                        offset = 0;
     //                                   start_flag = 0;
                                                        }
                                //fprintf(stderr,"offset is :%ld\n",offset);
			    }
			}else{continue;}//if (seq == 0 || start_flag == 1){
		}else{continue;} // if n>0;


        pthread_testcancel();
    } //while loop

    // Thread success!
    //fflush(stdout);
    return THREAD_OK;
}// static 

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
