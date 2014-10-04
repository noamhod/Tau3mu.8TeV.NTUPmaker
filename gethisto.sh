#!/bin/sh
if [ -z "$3" ] ; then
	echo "ERROR, you should provide [1] the root file name, [2] the histogram name, [3] the histogram type"
	echo "./gethisto.sh  /path/to/file.root  hisogramX  TH1"
	exit 0
fi
fname=$1
hname=$2
htype=$3
root -b -l -q gethisto.C++ --f=${fname} --h=${hname} --t=${htype}
