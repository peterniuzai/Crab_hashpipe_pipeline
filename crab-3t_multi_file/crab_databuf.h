#include <stdint.h>
#include <stdio.h>
#include "hashpipe.h"
#include "hashpipe_databuf.h"

#define CACHE_ALIGNMENT         64
#define N_INPUT_BLOCKS          3 
#define N_OUTPUT_BLOCKS         3

#define N_CHAN_PER_PACK		512		//number of channels per packet
#define N_PACKETS_PER_SPEC	8		//number of packets per spectrum
#define N_BYTES_DATA_POINT	2		// number of bytes per datapoint in packet
#define N_POLS_CHAN		4		//number of polarizations per channel
#define N_BYTES_COUNTER		8		// number bytes of counter
#define N_CHANS_SPEC		(N_CHAN_PER_PACK * N_PACKETS_PER_SPEC) // not including the poles in spec. if we have 4 polarition, the total chanels should time 4
#define DATA_SIZE_PACK		(unsigned long)(N_CHAN_PER_PACK * N_POLS_CHAN *  N_BYTES_DATA_POINT) //(4096 MBytes)doesn't include header for each packet size 
#define PKTSIZE			((N_CHAN_PER_PACK * N_POLS_CHAN *  N_BYTES_DATA_POINT) + N_BYTES_COUNTER)
#define N_BYTES_PER_SPEC	(DATA_SIZE_PACK*N_PACKETS_PER_SPEC)
#define N_SPEC_BUFF		512//128
#define BUFF_SIZE		(N_SPEC_BUFF*N_BYTES_PER_SPEC) // including all 4 polaration
#define N_CHANS_BUFF		(N_CHANS_SPEC*N_SPEC_BUFF*N_POLS_CHAN)
#define N_SPEC_PER_FILE		120320//960000//20096  // int(time(s)/T_samp(s)/N_SPEC_BUFF)*N_SPEC_BUFF   (20/0.001/128*128)
#define N_MBYTES_PER_FILE	(N_SPEC_PER_FILE * N_BYTES_PER_SPEC /1024/1024/ N_POLS_CHAN) //for now we only save I polaration into disk. Note Here We use Mbytes.
//extern unsigned long long miss_pkt;
// Used to pad after hashpipe_databuf_t to maintain cache alignment
typedef uint8_t hashpipe_databuf_cache_alignment[
  CACHE_ALIGNMENT - (sizeof(hashpipe_databuf_t)%CACHE_ALIGNMENT)
];

/* INPUT BUFFER STRUCTURES
  */
typedef struct crab_input_block_header {
   unsigned long long  mcnt;                    // mcount of first PACK
} crab_input_block_header_t;

typedef uint8_t crab_input_header_cache_alignment[
   CACHE_ALIGNMENT - (sizeof(crab_input_block_header_t)%CACHE_ALIGNMENT)
];

typedef struct crab_input_block {
   crab_input_block_header_t header;
   crab_input_header_cache_alignment padding; // Maintain cache alignment
   short int  data[N_CHANS_SPEC*N_SPEC_BUFF*N_POLS_CHAN];//*sizeof(char)];//512 FFT channels * 4 IFs * 2bytes = 4096Bytes
} crab_input_block_t;

typedef struct crab_input_databuf {
   hashpipe_databuf_t header;
   hashpipe_databuf_cache_alignment padding; // Maintain cache alignment
   crab_input_block_t block[N_INPUT_BLOCKS];
} crab_input_databuf_t;


/*
  * OUTPUT BUFFER STRUCTURES
  */
typedef struct crab_output_block_header {
   unsigned long long mcnt;
} crab_output_block_header_t;

typedef uint8_t crab_output_header_cache_alignment[
   CACHE_ALIGNMENT - (sizeof(crab_output_block_header_t)%CACHE_ALIGNMENT)
];

typedef struct crab_output_block {
   crab_output_block_header_t header;
   crab_output_header_cache_alignment padding; // Maintain cache alignment
//   uint64_t sum;
   short int I[N_CHANS_SPEC*N_SPEC_BUFF];
//   short int Q[N_CHANS_SPEC*N_SPEC_BUFF];
//   short int U[N_CHANS_SPEC*N_SPEC_BUFF];
//   short int V[N_CHANS_SPEC*N_SPEC_BUFF];

} crab_output_block_t;

typedef struct crab_output_databuf {
   hashpipe_databuf_t header;
   hashpipe_databuf_cache_alignment padding; // Maintain cache alignment
   crab_output_block_t block[N_OUTPUT_BLOCKS];
} crab_output_databuf_t;

/*
 * INPUT BUFFER FUNCTIONS
 */
hashpipe_databuf_t *crab_input_databuf_create(int instance_id, int databuf_id);

static inline crab_input_databuf_t *crab_input_databuf_attach(int instance_id, int databuf_id)
{
    return (crab_input_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int crab_input_databuf_detach(crab_input_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline void crab_input_databuf_clear(crab_input_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline int crab_input_databuf_block_status(crab_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_input_databuf_total_status(crab_input_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}

static inline int crab_input_databuf_wait_free(crab_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_input_databuf_busywait_free(crab_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_input_databuf_wait_filled(crab_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_input_databuf_busywait_filled(crab_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_input_databuf_set_free(crab_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_free((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_input_databuf_set_filled(crab_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_filled((hashpipe_databuf_t *)d, block_id);
}

/*
 * OUTPUT BUFFER FUNCTIONS
 */

hashpipe_databuf_t *crab_output_databuf_create(int instance_id, int databuf_id);

static inline void crab_output_databuf_clear(crab_output_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline crab_output_databuf_t *crab_output_databuf_attach(int instance_id, int databuf_id)
{
    return (crab_output_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int crab_output_databuf_detach(crab_output_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline int crab_output_databuf_block_status(crab_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_output_databuf_total_status(crab_output_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}

static inline int crab_output_databuf_wait_free(crab_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_output_databuf_busywait_free(crab_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_free((hashpipe_databuf_t *)d, block_id);
}
static inline int crab_output_databuf_wait_filled(crab_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_output_databuf_busywait_filled(crab_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_output_databuf_set_free(crab_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_free((hashpipe_databuf_t *)d, block_id);
}

static inline int crab_output_databuf_set_filled(crab_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_filled((hashpipe_databuf_t *)d, block_id);
}


