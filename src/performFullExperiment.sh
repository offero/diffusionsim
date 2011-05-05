#!/bin/bash

SIM_OUTPUT_PATH=/home/prima/Development/tmp/disim/100Trials31Nodes

python disim.py simulate --trials=100 --nodes=31 --direction=both -o $SIM_OUTPUT_PATH

