#!/usr/bin/env bash
source /etc/profile.d/modules.sh
module load xray
module load anaconda3/5.2
#set pathes
INPUTDATAFILE=$(realpath $1)
if [ $? -ne 0 ]; then
  exit 1
fi
INPUTDATAPATH="${INPUTDATAFILE%/*}"
OUTPUTDATAPATH="$INPUTDATAPATH/dozor"
if [[ $OUTPUTDATAPATH == "/asap3/petra3/gpfs/"* && $OUTPUTDATAPATH =~ "/raw/" ]]; then
  OUTPUTDATAPATH="/beamline/p11/current/raw/${OUTPUTDATAPATH##*/raw/}"
fi
if [[ $OUTPUTDATAPATH =~ "/raw/" ]]; then
  OUTPUTDATAPATH="${OUTPUTDATAPATH/\/raw\///processed/}"
fi
cd "$OUTPUTDATAPATH" || exit 1

#wait for master file
MAXATTEMPTS=30
ATTEMPT=0
h5stat "$INPUTDATAFILE" >/dev/null 2>&1
COMPLETE=$?
while [ $COMPLETE -ne 0 ]; do
  #echo $COMPLETE
  ATTEMPT=$[$ATTEMPT + 1]
  if [ "$ATTEMPT" -ge "$MAXATTEMPTS" ]; then
    echo "Cannot read master file: $INPUTDATAFILE"
    echo "$ATTEMPT failed attempts"
    exit 1
  fi
  sleep 1  # Use sleep because inotify-tools may not be available
  h5stat "$INPUTDATAFILE" >/dev/null 2>&1
  COMPLETE=$?
done

#retrieve meta data
SIZEX=$(h5dump -d /entry/instrument/detector/module/data_size "$INPUTDATAFILE" |awk '$1 ~ /\(0\)/ {print ($2+0);exit}')
SIZEY=$(h5dump -d /entry/instrument/detector/module/data_size "$INPUTDATAFILE" |awk '$1 ~ /\(0\)/ {print $3;exit}')
DETDIST=$(h5dump -d /entry/instrument/detector/detector_distance "$INPUTDATAFILE" |awk '$1 ~ /\(0\)/ {print 1000*$2;exit}'|bc)
WAVELENGTH=$(h5dump -d /entry/instrument/beam/incident_wavelength "$INPUTDATAFILE" |awk '$1 ~ /\(0\)/ {print $2;exit}')
ORGX=$(h5dump -d /entry/instrument/detector/beam_center_x "$INPUTDATAFILE" |awk '$1 ~ /\(0\)/ {print $2;exit}')
ORGY=$(h5dump -d /entry/instrument/detector/beam_center_y "$INPUTDATAFILE" |awk '$1 ~ /\(0\)/ {print $2;exit}')
if [ $SIZEX -gt "3000" ]; then
  IXMAX=2200
  IYMIN=2100
  IYMAX=2250
else
  IXMAX=1200
  IYMIN=1020
  IYMAX=1170
fi

echo " N    | SPOTS     Main     Visible"
echo "image | num.of   Score    Resolution"
echo "------------------------------------"

#loop through all container files and start dozor for each
CONTAINER=0
for DATACONTAINER in $(h5dump -n "$INPUTDATAFILE" | grep 'ext link' | awk '{print $5}') ; do
  MAXATTEMPTS=30
  ATTEMPT=0
  CONTAINER=$[$CONTAINER + 1]
  h5stat "$INPUTDATAPATH/$DATACONTAINER" >/dev/null 2>&1
  COMPLETE=$?
  while [ $COMPLETE -ne 0 ]; do
    #echo $COMPLETE
    ATTEMPT=$[$ATTEMPT + 1]
    if [ "$ATTEMPT" -ge "$MAXATTEMPTS" ]; then
      echo "Cannot read data file: $DATACONTAINER"
      echo "$ATTEMPT failed attempts"
      exit 1
    fi
    sleep 1  # Use sleep because inotify-tools may not be available
    h5stat "$INPUTDATAPATH/$DATACONTAINER" >/dev/null 2>&1
    COMPLETE=$?
  done

  #date +%s.%N
  START=$(h5dump -a /entry/data/data/image_nr_low "$INPUTDATAPATH/$DATACONTAINER" |awk '$1 ~ /\(0\)/ {print $2;exit}')
  STOP=$(h5dump -a /entry/data/data/image_nr_high "$INPUTDATAPATH/$DATACONTAINER" |awk '$1 ~ /\(0\)/ {print $2;exit}')
  NIMAGES=$[$STOP - $START +1]
  cat << EOF >"$OUTPUTDATAPATH/dozorr$CONTAINER.dat"
nx $SIZEX
ny $SIZEY
pixel 0.075
detector_distance $DETDIST
X-ray_wavelength $WAVELENGTH
orgx $ORGX
orgy $ORGY
ix_min 0
ix_max $IXMAX
iy_min $IYMIN
iy_max $IYMAX
first_image_number $START
number_images $NIMAGES
name_template_image $INPUTDATAFILE
library /opt/xray/XDSAPP/plugins/dectris-neggia.so
EOF
  sleep 1
  #date +%s.%N
  dozor "$OUTPUTDATAPATH/dozorr$CONTAINER.dat" |tee "$OUTPUTDATAPATH/dozorr$CONTAINER.out" |awk '($1 ~ /[[:digit:]]/) {print $0}' &
done

python3 Monitor_for_dozor_real_time_v1.py "$INPUTDATAPATH"

wait
echo "------------------------------------"

