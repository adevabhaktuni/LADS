#!/bin/bash
START=$(date +%s)
## Change to deisotope directory

INPUT=""
OUTPUT=""
LEVEL=2
DE_ONLY=true
LADS_ONLY=true
SCORE_ONLY=true
OUT_DEFINED=false
THREAD=""
PROFILE=false
DATABASE_SEARCH=""

while getopts "hd:o:l:r:tps:" opt;
do 
	case "$opt" in
		"h")
			echo "

	Deisotoper + LADS search

	Developed by Arun Devabhaktuni and the Elias Lab at Stanford University.

	Using this script will launch both deisitoper and LADS search. To run, 
	enter the following command:

	sh wrap.sh -d INPUT -o OUTPUT -l MS_LEVEL (optional) -r RESTRICTIONS (optional)

	where -d is the input mzXML file, -o is the desired return directory for LADS
	and LADS scoring, -l related so the MS level that deisotoper will search for,
	and -r is the single piece that you wish to run (optional). Input must be entered
	with full path. By default, if no output directory is given, this script will 
	default to storing all DTAs and LADS output files to the input directory. 
	If the MS level is not specified, then it will default to MS2. 

	Possible MS levels:
	1 = MS1
	2 = MS2
	3 = MS1 + MS2
	4 = MS3
	5 = MS1 + MS3
	6 = MS2 + MS3
	7 = MS1 + MS2 + MS3

	Possible restrictions with -r:
	deisotope = Only run deisotoper
	lads = Only run LADS search
	score = Only run LADS scoring function

	"
	exit;
	;;
	"d") 
		if [[ "$OPTARG" == /* ]]; then
			INPUT=$OPTARG
		else
			INPUT=$(readlink -f $OPTARG)
		fi
		echo "INPUT $INPUT"
	;;
	"o")
		if [[ "$OPTARG" == /* ]]; then
			OUTPUT=$OPTARG
		else
			OUTPUT=$(readlink -f $OPTARG)
		fi
		OUT_DEFINED=true
		echo "OUTPUT $OUTPUT"
	;;
	"l")
		LEVEL=$OPTARG
	;;
	"r")
		if [[ "$OPTARG" == deisotope ]]; then
			echo "Only executing deisotope function."
			LADS_ONLY=false
			SCORE_ONLY=false
		elif [[ "$OPTARG" == lads ]]; then
			echo "Only executing LADS function."
			DE_ONLY=false
			SCORE_ONLY=false
		elif [[ "$OPTARG" == score ]]; then
			echo "Only executing LADS scoring function."
			DE_ONLY=false
			LADS_ONLY=false
		fi
	;;
	"t")
		THREAD="Threaded"
	;;
	"p")
		PROFILE=true
	;;
	"s")
		DATABASE_SEARCH=$OPTARG
	;;
	esac
done

if ! $OUT_DEFINED; then
	echo "No output defined"
	IFS='/' read -a split <<< "$INPUT"
	full_file=${split[-1]}
	IFS='.' read -a file <<< "$full_file"
	file=${file[0]}
	OUTPUT=`pwd`/$file
	mkdir $OUTPUT
	echo "Output set $OUTPUT"
fi

if [[ ! -d $OUTPUT ]]; then
	mkdir $OUTPUT
fi

if $DE_ONLY; then
	if [[ "$OUTPUT" == */ ]]; then
		TMP=$OUTPUT'dta/'
	else
		TMP=$OUTPUT/'dta/'
	fi
elif [[ !$DE_ONLY ]]; then
	TMP=$INPUT'/'
elif [[ !$DE_ONLY ]] && [[ !$LADS_ONLY ]] && [[ $SCORE_ONLY ]]; then
	if [[ "$INPUT" == */ ]]; then
		TMP=$INPUT'dta/'
	else
		TMP=$INPUT/'dta/'
	fi
fi

mkdir $TMP
cd deisotope
echo ""

if $DE_ONLY ; then 
	echo "RUNNING DEISOTOPER ON PROVIDED MZXML"
	echo "LOADING FROM $INPUT" 
	echo "EXPORTING TO $TMP"
	php deisotoper_mzxml_conversion.php -n -f $INPUT -o $TMP -l $LEVEL
	DE=$(date +%s)
	DIFF=$(( $DE - $START ))
	echo "Deisotope took $DIFF seconds" 

fi

if [[ "$OUTPUT" == */ ]]; then
	OUTPUT=${OUTPUT%?}
fi

cd ../LADS/Scripts
if $LADS_ONLY ; then
	echo "RUNNING LADS ON RECENT DTAs"
	echo "LOADING FROM $TMP"
	pwd
	START=$(date +%s)

	if $PROFILE; then
		echo "Preparing to profile LADS"
		kernprof.py -l SequenceDTAsTDV$THREAD.py -p 5 -a 4 -A 12 -e 300 -E 1000 -l 0.9 -S ../Misc/symbolmap.txt -k ../Misc/TDVcols.txt -i ../Misc/LADS_LysC_guanStatic_lDimethDiff_hDimethDiff_cluster_light_heavy.ini -c ../Scoring_Functions/LysC_HCD_likelihood_prior_b_y_a_pos_AAClass_singleNeutralLoss_IntBin5_config.txt -m ../Scoring_Functions/ath013833_ath013837_LysC_guanDiMeth_LogTIC_intBin5_bya_singNeutLoss_AAClass_pos_likelihood_prior.model -d $TMP -o $OUTPUT/output.txt | tee $OUTPUT/output.log
	else 

	python2.6 SequenceDTAsTDV$THREAD.py -p 5 -a 4 -A 12 -e 300 -E 1000 -l 0.9 -S ../Misc/symbolmap.txt -k ../Misc/TDVcols.txt -i ../Misc/LADS_LysC_guanStatic_lDimethDiff_hDimethDiff_cluster_light_heavy.ini -c ../Scoring_Functions/LysC_HCD_likelihood_prior_b_y_a_pos_AAClass_singleNeutralLoss_IntBin5_config.txt -m ../Scoring_Functions/ath013833_ath013837_LysC_guanDiMeth_LogTIC_intBin5_bya_singNeutLoss_AAClass_pos_likelihood_prior.model -d $TMP -o $OUTPUT/output.txt | tee $OUTPUT/output.log
	
	fi
	
	LADS=$(date +%s)
	DIFF=$(( $LADS - $START ))
	echo "LADS search took $DIFF seconds" 
fi

if $SCORE_ONLY ; then
	echo "RUNNING LADS SCORING FUNCTION"
	pwd
	echo "LOADING FROM $OUTPUT with $TMP as DTA source"

	if [[ "$OUTPUT" != "" ]]; then
		$DATABASE_SEARCH = " -s $DATABASE_SEARCH"
	fi

	START=$(date +%s)
	python2.6 PostLADResultsAndWriteAccuracySVM.py -i ../Misc/LADS_LysC_guanStatic_lDimethDiff_hDimethDiff_cluster_light_heavy.ini -c ../Scoring_Functions/LysC_HCD_likelihood_prior_b_y_a_pos_AAClass_singleNeutralLoss_IntBin5_config.txt -m ../Scoring_Functions/ath013833_ath013837_LysC_guanDiMeth_LogTIC_intBin5_bya_singNeutLoss_AAClass_pos_likelihood_prior.model -S ../Misc/symbolmap.txt -D ../Scoring_Functions/20130204_94Feature_SVMFormat_14Datasets_VB_guanDimeth_Tryp_LysC_removeNAN.model -L $OUTPUT/output.log -o $OUTPUT/output_PostScore.tdv -d $TMP $DATABASE_SEARCH

	SCORE=$(date +%s)
	DIFF=$(( $SCORE - $START ))
	echo "Scoring took $DIFF seconds" 
fi
TOT=$(date +%s)
DIFF=$(( $TOT - $START ))
echo "Total run time $DIFF"
echo "OPERATIONS COMPLETE"

exit
