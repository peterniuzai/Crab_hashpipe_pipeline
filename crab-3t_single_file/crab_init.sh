#!/bin/bash
#./send_packet.py &
hashpipe -p crab_hashpipe -I 0 -c 1 crab_net_thread -c 2 crab_gpu_thread -c 3 crab_output_thread

