#!/bin/bash
#=============================================================================
INITALS="GRC";
MAX=22;
WGETOP=" --trust-server-names -q ";
ONWAY="ftp://130.14.250.12/genomes/Homo_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38.p2_chr"
#-----------------------------------------------------------------------------
function downloadEach
  {
  PATTERN="unexpected";
  for((;;)); 
    do
    sleep 2;
    xtmp="`wget $1 $2 -O $4-X$3`"
    if [[ "$xtmp" == *"$PATTERN"* ]];  
      then
      echo "Unexpected EOF found, re-downloading C$3 ...";
      continue;
    else
      echo "wget stderr:$xtmp";
      echo "Downloaded $4 C$3 with success!";
      break;
    fi
    done
  }

echo "Downloading and filtering $INITALS sequences ..."
for((x=1 ; x <= $MAX ; ++x)); 
  do
  ZPATH="$ONWAY$x.fa.gz";
  downloadEach "$WGETOP" "$ZPATH" "$x" "$INITALS";
  zcat $INITALS-X$x > $INITALS$x;
  echo "$INITALS C$x filtered!";
  done

CHR=23;
FIELD="X";
ZPATH="$ONWAY$FIELD.fa.gz";
downloadEach "$WGETOP" "$ZPATH" "$CHR" "$INITALS";
zcat $INITALS-X$CHR > $INITALS$CHR;
echo "$INITALS CX filtered";

CHR=24;
FIELD="Y";
ZPATH="$ONWAY$FIELD.fa.gz";
downloadEach "$WGETOP" "$ZPATH" "$CHR" "$INITALS";
zcat $INITALS-X$CHR > $INITALS$CHR;
echo "$INITALS CY filtered";

cat HS1 HS2 HS3 HS4 HS5 HS6 HS7 HS8 HS9 HS10 HS11 HS12 HS13 HS14 HS15 HS16 /
HS17 HS18 HS19 HS20 HS21 HS22 HS23 HS24 > HS;
rm *GRC-* -f
echo "Done!"
#=============================================================================
