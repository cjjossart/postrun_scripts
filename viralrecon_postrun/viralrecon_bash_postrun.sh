##!/bin/bash

echo "-- Starting post-run bash script --"
echo "Working directory: $PWD"

echo "Copying Nextflow out log file..."
cp ${NXF_OUT_FILE} log.txt

outdir_cmd="export outdir=`cat log.txt | grep "outdir" | awk '{print $NF}'`"
echo "Getting 'outdir' from log file and exporting as a variable with command: '$outdir_cmd'"
$outdir_cmd
echo "Variable 'outdir' is: $outdir"

echo "Installing pysam dependencies"
echo "Command: yum install autoconf automake make gcc perl-Data-Dumper zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel ncurses-devel python3-devel"
yum -y install autoconf automake make gcc perl-Data-Dumper zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel ncurses-devel python3-devel
echo "Installing pandas boto3 cython and pysam"
pip3 install pandas boto3 cython pysam

cmd="python3 viralrecon_postrun.py $outdir"
echo "-- Running python script with command '$cmd'... --"
echo " "
$cmd
