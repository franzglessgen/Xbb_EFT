#!/bin/bash

STARTTIME=$(date +%s.%N)
ulimit -c 0
ulimit -S -c 0

export TERM=""

echo "python dummy.py"
eval "python ./dummy.py"
STARTTIME=$(date +%s.%N)
STARTTIME=$(date +%s.%N)

EXITCODE=$?
echo "--------------------------------------------------------------------------------"
echo "exit code: $EXITCODE"
echo "EndTime: "`date +"%Y-%m-%d %T"`
ENDTIME=$(date +%s.%N)
DIFFTIME=$(echo "($ENDTIME - $STARTTIME)/60" | bc)
echo "duration (real time): $DIFFTIME minutes"
if [ "$EXITCODE" -ne "0" ]; then
    if [ "$noretry" = "1" ]; then
        echo "--- STOP ---"
    else
        echo "--- RETRY ---"
    fi
else
    echo "--- OK ---"
fi
echo
echo "Exiting runAll.sh"
echo
exit $EXITCODE
