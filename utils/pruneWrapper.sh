#! /bin/bash

PRUNE="../rnr-prune"

NUM=$1
RNR_FILE=$2
BS_FILE=$3
TRE_FILE=$4


if [ $# -lt 3 ]; then
    echo -e  "This is a tiny helper script for obtaining pruned trees from a RogueNaRok search. \n
Call it as follows:\n
./pruneWrapper.sh <numToPrune> <RogueNaRok_droppedRogues-file> <bootstrap-tree-file>  <ML-tree-file>\n
* <numToPrune> prune the first x rogues that RogueNaRok identified (the order in the output is important). \n  If you want to prune all rogues as identified by RogueNaRok, simply enter an arbitrarily large number in this position \n
* further arguments are self-explainatory.\n
* output is written to stdout\n" 
    exit
fi

id=$(echo $RNR_FILE  |awk -F . '{print $NF}')

excludeFile=$id.$NUM.exclude

cat $RNR_FILE | tail -n +3  | head -n $NUM | cut -f 3  > $excludeFile



if [ "$TRE_FILE" ==  "" ]; then
    $PRUNE -x $excludeFile -i  $BS_FILE  -n $id > /dev/null
    if [ $? == 0 ]; then
	echo "wrote file RnR-prune_prunedBootstraps.$id"
    else
	echo "something went wrong, maybe file RnR-prune_prunedBootstraps.$id already exists (delete it to re-run)"
    fi
    
else
    $PRUNE -x $excludeFile -i  $BS_FILE -t $TRE_FILE  -n $id  > /dev/null
    if [ $? == 0 ]; then
	echo "wrote file RnR-prune_prunedBootstraps.$id and RnR-prune_prunedBestTree.$id"
    else
	echo "something went wrong, maybe file RnR-prune_prunedBootstraps.$id or RnR-prune_prunedBestTree.$id already exist (delete them to re-run)"
    fi
    
fi







