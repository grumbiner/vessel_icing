#!/bin/sh
set -e

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
source=https://ftpprd.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs.${tag}/${cyc}

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
    
#set -xe

hr=000 
#Get static fields:
   #sst
  if [ ! -f sst.$res.$tag.f$hr.grib2 ] ; then
     ${BIN_DIR}/get_inv.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}.idx | grep -f ${BIN_DIR}/sst_parms | \
    ${BIN_DIR}/get_grib.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}        sst.$res.$tag.f$hr.grib2
  fi

   #land, ice
  if [ ! -f landice.$res.$tag.f$hr.grib2 ] ; then
   ${BIN_DIR}/get_inv.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}.idx | grep -f ${BIN_DIR}/landice_parms | \
  ${BIN_DIR}/get_grib.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}        landice.$res.$tag.f$hr.grib2
  fi

#240 hours is the limit of 3-hourly output from GFS
while [ $hr -le 240 ]
do
   # running inputs -- u10, v10, t2m
  if [ ! -f running.$res.$tag.f$hr.$cyc.grib2 ] ; then
   ${BIN_DIR}/get_inv.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}.idx | grep -f ${BIN_DIR}/running_parms | \
  ${BIN_DIR}/get_grib.pl $source/gfs.t${cyc}z.pgrb2.${res}.f${hr}        running.$res.$tag.f$hr.$cyc.grib2
  fi

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
if [ ! -f sst ] ; then
  ${BIN_DIR}/wgrib2 sst.$res.$tag.f$hr.grib2 | ${BIN_DIR}/wgrib2 -i sst.$res.$tag.f$hr.grib2 -no_header -order we:ns -bin sst
fi
if [ ! -f landice ] ; then
  ${BIN_DIR}/wgrib2 landice.$res.$tag.f$hr.grib2 | ${BIN_DIR}/wgrib2 -i landice.$res.$tag.f$hr.grib2 -no_header -order we:ns -bin landice
fi

if [ ! -f running_input ] ; then
  cat running.*.grib2 > all.grib2
  ${BIN_DIR}/wgrib2 all.grib2 | ${BIN_DIR}/wgrib2 -i all.grib2 -no_header -order we:ns -bin running_input
  #Preliminary cleanup 
  rm running.*.grib2
fi

########## # Run the model ################################################################
cp $MODEL_DIR/*.py .
#time python3 -m cProfile -o prof.stat icing.py > icing.out
#python3 statview.py > stats.out

time python3 icing.py > icing.out


######### # Make some output #############################################################
# -- now in python, via matplotlib

#grep hist icing.out > hist
#cp $MODEL_DIR/gnuin .
#gnuplot < gnuin
mv *.png ~/grumbinescience.org/icing/

########## # Make some output #############################################################
#Final cleanup
rm all.grib2 running_input
