#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
cmsswDir="/uscmst1b_scratch/lpc1/lpcmuon/zhenhu/FNAL/dimuon/CMSSW_7_4_6_patch6/src/"
inputFiles=""
dcacheDirIn="/eos/uscms/store/user/zhenhu/MuOnia/Onia2MuMuRootuple-Run2015B-MuOina-v3/f64ae4aafa8965110120322f8d8de3c5"

j=0
for files in `ls $dcacheDirIn | grep Rootuple_`
do 
	inputFiles=$files;
	echo $inputFiles;
	jobNb=${j};
	let j=${j}+1;
	command="runRemoveDuplicate_${jobNb}_C.so,runRemoveDuplicate_${jobNb}_C.d,runRemoveDuplicate_${jobNb}_C_ACLiC_dict_rdict.pcm";
	jobCScript="runRemoveDuplicate_${jobNb}.C";
	scriptName="Run_${jobNb}.csh";
	condorScriptName="runOnCondor_${jobNb}";
	#cat myntuple.h | sed "s/NUMBER/${jobNb}/g" > ${anaHeader};
	#cat myntuple.C | sed "s/NUMBER/${jobNb}/g" > ${anaCScript};
	cat runRemoveDuplicate.C | sed "s_INPUTPATH_${dcacheDirIn}_" | sed "s-INPUTFILE-${inputFiles}-" | sed "s/NUMBER/${jobNb}/g" > ${jobCScript};
	cat Run.csh | sed "s-FILENAME-${jobCScript}-" | sed "s-CMSSWDIR-${cmsswDir}-" > ${scriptName};
	chmod +x ${scriptName}
	cat runOnCondor | sed "s/SCRIPT/${scriptName}/" | sed "s-COMMAND-${command}-" | sed "s/JOBC/${jobCScript}/" > ${condorScriptName}
	
	root -l <<EOF
.L ${jobCScript}++
.q
EOF
	
	condor_submit ${condorScriptName}
done
