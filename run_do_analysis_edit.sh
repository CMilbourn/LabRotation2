# !/bin/bash
# run_do_analysis.sh sourcedata 

#Prints out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 <sourcedata folder> e.g. ./run_do_analysis_edit.sh /Users/colette/sourcedata/"; 
 #exit 1 ;
fi

# Input one is the sourcedata directory path
src=$1 

for n in {3..9}
do
	subj_id=sub-`printf "%02d\n" $n`
	echo "Subject ID:"$subj_id
	
	echo "Running do_analysis_edit.sh: ${src} ${subj_id}"
	${codedir}./do_analysis_edit.sh ${src} ${subj_id} 
done