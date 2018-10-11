#!/bin/sh
#set -x

#Get data and run vessel icing model
#Robert Grumbine 11 October 2018

#General notes:
#  As currently structured, a model run requires approximately 
#    1.5 Gb running space
#    3:30 to download the approx 177 Mb of grib2 data (at home)
#    7:30 to run the python code.


export tag=${tag:-`date +"%Y%m%d"`}
export cyc=${cyc:-00}

#res=0p50, 1p00 (0.25, 0.50, 1.00 degree lat-long grids from the NWS GFS model)
export res=${res:-0p25}

# SOURCE
#The sources provide inputs through curl calls invoked by get_grib.pl and get_inv.pl
#ftpprd is NWS operational (and shorter retention)
#nomads is developmental, but longer retention
#cfsv2 is the climate forecast system reanalysis, provinding information 1979-present
source=http://ftpprd.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs.${tag}${cyc}
#source=http://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.${tag}${cyc}
#cfsv2 a:
#cfsv2 b:

BASE=`pwd`
BIN_DIR=${BIN_DIR:-${BASE}/../bin}
if [ ! -d $BIN_DIR ] ; then
  echo no BIN_DIR at $BIN_DIR exiting
  exit
fi

MODEL_DIR=${MODEL_DIR:-$BASE}
if [ ! -d $MODEL_DIR ] ; then
  echo no MODEL_DIR at $MODEL_DIR exiting
  exit
fi

#cd to a working directory
RUN_DIR=${RUN_DIR:-running}
if [ ! -d $RUN_DIR ] ; then
  mkdir -p $RUN_DIR
  err=$?
  if [ err -ne 0 ] ; then
    echo Making rundir failed with error $err exiting now
    exit
  fi
fi
cd $RUN_DIR
    

hr=000 
#Get static fields:
   #sst
   ${BIN_DIR}/get_inv.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}.idx | grep -f ${BIN_DIR}/sst_parms | \
  ${BIN_DIR}/get_grib.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}        sst.$res.$tag.f$hr.grib2

   #land, ice
   ${BIN_DIR}/get_inv.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}.idx | grep -f ${BIN_DIR}/landice_parms | \
  ${BIN_DIR}/get_grib.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}        landice.$res.$tag.f$hr.grib2

#240 hours is the limit of 3-hourly output from GFS
while [ $hr -le 240 ]
do
   # running inputs -- u10, v10, t2m
   ${BIN_DIR}/get_inv.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}.idx | grep -f ${BIN_DIR}/running_parms | \
  ${BIN_DIR}/get_grib.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}        running.$res.$tag.f$hr.$cyc.grib2

  hr=`expr $hr + 3`
  if [ $hr -le 10 ] ; then
    hr=0$hr
  fi
  if [ $hr -le 100 ] ; then
    hr=0$hr
  fi

done

########## Generate binaries for model (once pygrib is working properly, this can be skipped)
hr=000
wgrib2 sst.$res.$tag.f$hr.grib2 | wgrib2 -i sst.$res.$tag.f$hr.grib2 -no_header -order we:ns -bin sst
wgrib2 landice.$res.$tag.f$hr.grib2 | wgrib2 -i landice.$res.$tag.f$hr.grib2 -no_header -order we:ns -bin landice
cat running.*.grib2 > all.grib2
wgrib2 all.grib2 | wgrib2 -i all.grib2 -no_header -order we:ns -bin running_input

########## # Run the model ##################################################################
cp $MODEL_DIR/*.py .
time python icing.py > icing.out
