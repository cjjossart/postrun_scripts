##!/bin/bash

echo "-- Starting post-run script --"
echo "Working directory: $PWD"

echo "Copying Nextflow out log file..."
cp ${NXF_OUT_FILE} log.txt

outdir_cmd="export outdir=`cat log.txt | grep "outdir" | awk '{print $NF}'`"
echo "Getting 'outdir' from log file and exporting as a variable with command: '$outdir_cmd'"
$outdir_cmd
echo "Variable 'outdir' is: $outdir"

echo "Copying viralrecon wrapper script to current dir.."
aws s3 cp s3://dev-wslh-sequencing-analyses/scripts/viralrecon_wrapper.py .

echo "Installing pandas boto3 and cython"
echo " "
pip3 install pandas boto3 cython

echo " "
echo "Installing pysam dependencies"
echo "Command: yum install autoconf automake make gcc perl-Data-Dumper zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel ncurses-devel python3-devel"
yes | yum install autoconf automake make gcc perl-Data-Dumper zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel ncurses-devel python3-devel
echo "Installing pysam"
pip3 install pysam

cmd="python3 viralrecon_wrapper.py $outdir"
echo "-- Running python script with command '$cmd'... --"
echo " "
$cmd
