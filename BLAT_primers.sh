#!/bin/bash
#BLAT_primers.sh

echo "Hi, $USER. Compare primers for off-site hybridization."
echo

#BLAT required before proceeding

#convert all necessary files to 2bit format

DATADIR=db_oligo
AMPPRIMERDIR=ampprimers-skpp15
FILES="$DATADIR/*.fasta"
for f in $FILES
do
  #echo "Processing $f file..."
  result_string="${f/fasta/2bit}"
  ./faToTwoBit $f $result_string
done

FILES="$DATADIR/*.genes"
for f in $FILES
do
  #echo "Processing $f file..."
  result_string="${f/genes/2bit}"
  ./faToTwoBit $f $result_string
done

#screen amplification primers at oligo level
./gfServer -tileSize=6 -stepSize=1 -minMatch=1 -maxGap=3 -canStop start localhost 17779 db_oligo/LibAll_finaloligos_noAmp.2bit > /dev/null 2>&1 &
sleep 8
./gfClient -minScore=0 -minIdentity=0 -nohead localhost 17779 ./ $AMPPRIMERDIR/ampAll_200_selected.fasta DHFR_psl/ampAll_200_hits.psl
./gfServer stop localhost 17779
sleep 8
./pslPretty DHFR_psl/ampAll_200_hits.psl $DATADIR/LibAll_finaloligos_noAmp.2bit $AMPPRIMERDIR/ampAll_200_selected.fasta DHFR_psl/ampAll_200_hits_pretty.out

declare -a libnum=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15")
declare -a libnum2=("16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30")
declare -a oligonum=("4" "4" "4" "4" "4" "4" "4" "4" "4" "4" "4" "4" "4" "5" "5")

declare -a asmprim=("504" "504" "504" "504" "504" "504" "504" "504" "504" "504" "504" "504" "504" "504" "504")
declare -a asmprim2=("504" "504" "504" "504" "504" "504" "504" "504" "504" "504" "504" "504" "504" "504" "504")

COUNTER=0
## now loop through each other sample and BLAT each Assembly Primer pair against its Lib gene level
for i in "${libnum[@]}"
do
    
    ./gfServer -tileSize=6 -stepSize=1 -minMatch=1 -maxGap=3 -canStop start localhost 17779 $DATADIR/DHFR_Lib"${libnum[$COUNTER]}"_"${oligonum[$COUNTER]}"oligo.full_nRE_nPrim.2bit > /dev/null 2>&1 &
    sleep 4
    ./gfClient -minScore=0 -minIdentity=0 -nohead localhost 17779 ./ asmprimers-skpp20/skpp-"${asmprim[$COUNTER]}".fasta DHFR_psl/asm"${libnum[$COUNTER]}".psl
    ./gfServer stop localhost 17779
    sleep 4
    ./gfServer -tileSize=6 -stepSize=1 -minMatch=1 -maxGap=3 -canStop start localhost 17779 $DATADIR/DHFR_Lib"${libnum[$COUNTER]}"_"${oligonum[$COUNTER]}"oligo.codon2.full_nRE_nPrim.2bit > /dev/null 2>&1 &
    sleep 4
    ./gfClient -minScore=0 -minIdentity=0 -nohead localhost 17779 ./ asmprimers-skpp20/skpp-"${asmprim2[$COUNTER]}".fasta DHFR_psl/asm"${libnum2[$COUNTER]}".psl
    ./gfServer stop localhost 17779
    sleep 4
    ./pslPretty DHFR_psl/asm"${libnum[$COUNTER]}".psl $DATADIR/DHFR_Lib"${libnum[$COUNTER]}"_"${oligonum[$COUNTER]}"oligo.full_nRE_nPrim.2bit asmprimers-skpp20/skpp-"${asmprim[$COUNTER]}".fasta DHFR_psl/asm"${libnum[$COUNTER]}"_pretty.out
    ./pslPretty DHFR_psl/asm"${libnum2[$COUNTER]}".psl $DATADIR/DHFR_Lib"${libnum[$COUNTER]}"_"${oligonum[$COUNTER]}"oligo.codon2.full_nRE_nPrim.2bit asmprimers-skpp20/skpp-"${asmprim2[$COUNTER]}".fasta DHFR_psl/asm"${libnum2[$COUNTER]}"_pretty.out
    
    let COUNTER=COUNTER+1 
done