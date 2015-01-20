#!/bin/sh
if [ -z "$8" ] ; then
	echo "ERROR, you should provide the following arguments:"
	echo "  [1] build mode"
	echo "  [2] proc type"
	echo "  [3] run type"
	echo "  [4] output directory"
	echo "  [5] channel name"
	echo "  [6] master triplet type"
	echo "  [7] selection method"
	echo "  [8] split to N cores (OPTIONAL)"
	echo "e.g.:"
	echo "./execute.sh  build/buildrun/loadrun \ "
	echo "              NTUPmaker / [...] \ "
	echo "              analysis / skim \ "
	echo "              \$HOME/data/bout/ \ "
	echo "              period[A,B,C,D,E,G,H,I,J,L/all/[+MC]]/bbTomu15/Wtaunu_3mu... \ "
	echo "              muons / muid \ "
	echo "              CUTS:ALL / MVA:BDTG \ "
	echo "              N:j[N>=1, 0<j<=N] where N=number of jobs to split into and j is one of these subjobs"
	exit 0
fi

mode=$1
type=$2
run=$3
outdir=$4
chnl=$5
master=$6
method=$7
split=$8
if [[ "$1" == *build* ]] ; then
	echo "===> Note: I'm cleaning all libs !!!";
	./clean.sh
fi
root -b -l -q  compile.C  --mode=${mode}  --type=${type}  --run=${run}  --outDir=${outdir}  --chnl=${chnl}  --master=${master}  --method=${method}  --split=${split}