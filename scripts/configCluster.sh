#! /bin/bash
set -o nounset

#
# options
#

OUTFILE=submit_job.sh

# EXECMACHINE=localhost
EXECMACHINE=isdc

[ "$EXECMACHINE" == "isdc" ] && EXECHOME=/home/isdc/coupon # for remote machines
[ "$EXECMACHINE" == "localhost" ] && EXECHOME=$HOME # for remote machines

function writeJobFile
{

  JOBNAME=$1
  MEM=8G

  # build SGE script
  echo "#!/bin/bash" > $OUTFILE
  echo "#$ -N $JOBNAME" >> $OUTFILE
  echo "##$ -m bes" >> $OUTFILE
  echo "##$ -M jean.coupon@unige.ch" >> $OUTFILE
  echo "#$ -V" >> $OUTFILE
  echo "#$ -q test" >> $OUTFILE
  echo "##$ -pe ompi 10" >> $OUTFILE
  echo "#$ -cwd">> $OUTFILE
  echo "#$ -o /home/isdc/coupon/data/tmp/SGE_logs" >> $OUTFILE
  echo "#$ -e /home/isdc/coupon/data/tmp/SGE_logs" >> $OUTFILE
  echo "#$ -l h_rt=48:00:00 " >> $OUTFILE
  echo "#$ -l h_vmem=$MEM" >> $OUTFILE
  echo "#$ -l h=!isdc-cn16.astro.unige.ch" >>  $OUTFILE

  # echo "unset module" >> $OUTFILE

  local CMD="
    echo \$HOSTNAME
    date

    if [ \"\$#\" -ne 0 ]; then
      SGE_TASK_ID=\$1
    fi "

  echo "$CMD" >> $OUTFILE

  echo $OUTFILE

  return
}

function submitJob
{

  EXECDIR=tmp

  # convert paths
  sed -i ''  "s|$HOME|$EXECHOME|g" $JOBFILE

  echo "submitting $COUNTER jobs..."

  #
  # upload script ...
  #

  ssh $EXECMACHINE "mkdir -p $EXECDIR"
  scp $JOBFILE $EXECMACHINE:$EXECDIR

  scp scripts/photoz.py $EXECMACHINE:$EXECDIR

  #
  # ... and run it
  #

  # [ "$EXECMACHINE" = "canfar" ] && ssh $EXECMACHINE "chmod +x $EXECDIR/$(basename $EXEC); canfar_submit $EXECDIR/$(basename $JDL) ubuntu-server-14.04_hscPipe_Dec24_2016 c1-7.5gb-30; exit"
  # [ "$EXECMACHINE" == "picholine" ] && ssh $EXECMACHINE "cd $EXECDIR; for((i=0;i<$COUNTER;i++)); do bash $(basename $EXEC) \$i; done; exit;"
  [ "$EXECMACHINE" == "isdc" ] && ssh $EXECMACHINE "cd $EXECDIR; qsub -t 1-$COUNTER $(basename $JOBFILE); exit"
  [ "$EXECMACHINE" == "localhost" ] && ( cd $EXECHOME/$EXECDIR; for((i=1;i<=$COUNTER;i++)); do bash $(basename $JOBFILE) $i; done; )

  return
}
