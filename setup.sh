#!/bin/bash

. /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6/setup.sh
. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.13/x86_64-slc6-gcc47-opt/root/bin/thisroot.sh

export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/bin:$PATH
export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/lib:$LD_LIBRARY_PATH


# Now start it by typing root.
# To get 32bit binaries, simply replace x86_64 by i686 in above directories.
