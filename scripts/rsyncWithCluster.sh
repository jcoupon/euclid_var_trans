#! /bin/bash
set -o nounset

#
# synchronise with cluster ($1)
#

# ---------------------------------------------- #
# Help message
# ---------------------------------------------- #
if [ $# -lt 2 ]; then

  printf "\n\n     Synchronises with clusters\n\n"
  printf "Usage: $0 cluster [data,data_back]\n"
  exit
fi


if [ "$2" == "data" ]; then

  HOST_REMOTE=$1

  CURRENT=$HOME/data/euclid/varTrans
  CURRENT_REMOTE=data/euclid/varTrans
  echo $CURRENT_REMOTE

  # create the directory if it doesn't exist
  ssh $HOST_REMOTE "mkdir -p $CURRENT_REMOTE; exit"

  # synchronise with local directory
  # (don't follow links and don't copy the source, simply copy the link)
  echo "running rsync..."
  rsync -rlvzcup  --exclude 'OLD*' $CURRENT/ $HOST_REMOTE:$CURRENT_REMOTE
  # --delete

fi

if [ "$2" == "data_back" ]; then

  HOST_REMOTE=$1

  CURRENT=$HOME/data/euclid/varTrans
  CURRENT_REMOTE=data/euclid/varTrans
  echo $CURRENT_REMOTE

  # create the directory if it doesn't exist
  ssh $HOST_REMOTE "mkdir -p $CURRENT_REMOTE; exit"

  # synchronise with local directory
  # (don't follow links and don't copy the source, simply copy the link)
  echo "running rsync..."
  rsync -rlvzcup  --exclude 'OLD*' $HOST_REMOTE:$CURRENT_REMOTE/ $CURRENT
  # --delete

fi
