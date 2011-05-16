#!/bin/bash

SIM_OUTPUT_PATH=/home/prima/Development/tmp/disim/100Trials21Nodes-take2

export PYTHONPATH=$PYTHONPATH:`pwd`/disim
python disim/disim.py simulate --trials=100 --nodes=21 --direction=both -o $SIM_OUTPUT_PATH

