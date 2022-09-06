#!/usr/bin/env bash
cd /home/noah/FEMTA/spacebound/blue-origin/

rm -f auto.log
tmux new-session -d -s origin 'sudo chrt -f 16 ./origin.x file=experiments/freeze.e | tee auto.log'
