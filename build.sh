#!/bin/bash
status=0

STARTINGDIR=${PWD}
echo $STARTINGDIR
# cd cogFiles
PYTHONPATH=${PWD}/cogFiles
echo "WARNING: Replacing EmJetHistos.cc"
$HOME/.local/bin/cog.py -r ${STARTINGDIR}/EmJetHistos.cc
status=$?
if [ $status -ne 0 ] # If previous command was not successful
then
    cd $STARTINGDIR
    exit $status
fi
make clean
make
cd $STARTINGDIR
# exit $status
