#!/usr/bin/env python3

import boto3
import pandas as pd
import sys
import pysam
from os.path import exists


outdir = sys.argv[1]
outbucket = outdir.strip().replace("s3://", "").split("/")[0].strip()
prefix = "/".join(outdir.replace("s3://", "").split(outbucket)[-1].split("/")[1:])

if not prefix.endswith("/"):
    prefix = prefix + "/"

s3 = boto3.client('s3')

def get_summary_file(outdir, outfile):
    objects = s3.list_objects(Bucket=outbucket, Prefix=prefix)
    for object in objects["Contents"]:
        file_object = object['Key'].find("summary_variants_metrics_mqc.csv")
        if file_object < 0:
            pass
        else:
            obj = s3.get_object(Bucket=outbucket, Key=object["Key"])
            obj_df = pd.read_csv(obj["Body"])
            obj_df.to_csv(outfile, sep=',', index=False)
    return outfile


def get_bam_files(outdir):
    # get all bam files in the run for running samtools downstream
    bam_list = []
    objects = s3.list_objects(Bucket=outbucket, Prefix=prefix)
    for object in objects["Contents"]:
        ivar_bam = object['Key'].find("ivar_trim")
        if ivar_bam < 0:
            pass
        else:
            if object['Key'].endswith(".bam"):
                key = object['Key'].split("/")[-1]
                s3.download_file(outbucket, object['Key'], key)
                bam_list.append(key)
    return bam_list


def get_pangolin_data(int_outfile, fin_outfile):
    # takes all of the pangolin.csv files into a list and puts them into a combined dataframe
    df_list = []
    objects = s3.list_objects(Bucket=outbucket, Prefix=prefix)
    for object in objects["Contents"]:
        pangolin=object['Key'].find("pangolin")
        if pangolin < 0:
            pass
        else:
            if object['Key'].endswith(".csv"):
                obj = s3.get_object(Bucket=outbucket, Key=object["Key"])
                obj_df = pd.read_csv(obj["Body"])
                df_list.append(obj_df)
                pangolin_df = pd.concat(df_list)
                pangolin_df.to_csv(int_outfile, index=False)
    
    """ in order to generate the final postrun output, we need to be able to join dataframes on sample name. To do this, we need to
    parse out the reference fasta file name from the sample name of the pangolin output.
    """

    with open(int_outfile, 'r') as io, open(fin_outfile, 'w') as fo:
        line_dict = {}
        for line in io:
            if line.startswith("taxon"):  # keep the header line for our final combined pangolin output file
                header = line
            else:
                p_sample = line.split(",")[0]
                sample = p_sample.split(" ")[0]
                line_dict[sample] = ",".join(line.split(",")[1:])  # the dict has desired sample ID as the key and rest of line as value
        fo.write(header)
        for k, v in line_dict.items():
            fo.write(k + "," + v)
    return fin_outfile


def samtools_coverage(out_report, bam_list):
    cov_file_list = []
    cov_df_list = []
    # run samtools coverage on each BAM file. This gets "depth after trimming" and "1X coverage after trimming"
    for bam in bam_list:
        samplename = bam.split(".bam")[0].rstrip()
        coverage_file = samplename + ".cov.txt"
        try:
            print(f"Running samtools coverage on {samplename}...")
            pysam.coverage(bam, "-o", coverage_file)
            print(f"Samtools coverage for {samplename} successful.." + "\n")
            cov_file_list.append(coverage_file)
        except:
                print(f"Samtools coverage for {samplename} unsuccessful.." + "\n")

    print("Concatenating coverage reports...")
    try:
        for c_file in cov_file_list:
                sample_id = c_file.split(".")[0]
                cfile_df = pd.read_csv(c_file, sep="\t")
                cfile_df.insert(0, "Sample", sample_id)
                cov_df_list.append(cfile_df)
                coverage_df = pd.concat(cov_df_list)
        coverage_df.to_csv(out_report, index=False)
        print(f"Successfully combined coverage reports. Report: {out_report}")
    except:
        print(f"Combining coverage reports was unsuccessful.")
    return out_report
        

def get_nextclade_data(int_file, outfile):
    # same strategy as 'pangolin' above
    nextclade_list = []
    objects = s3.list_objects(Bucket=outbucket, Prefix=prefix)
    for object in objects["Contents"]:
        file_object = object['Key'].find("nextclade")
        if file_object < 0:
            pass
        else:
            if object['Key'].endswith(".csv"):
                obj = s3.get_object(Bucket=outbucket, Key=object["Key"])
                obj_df = pd.read_csv(obj["Body"], sep=";")
                nextclade_list.append(obj_df)
                nextclade_df = pd.concat(nextclade_list)
                nextclade_df.to_csv(int_file, index=False)
    with open(int_file, 'r') as io, open(outfile, 'w') as fo:
        line_dict = {}
        for line in io:
            if line.startswith("seqName"):
                header = line
            else:
                sample = line.split(" ")[0]
                line_dict[sample] = ",".join(line.split(",")[1:]) 
        fo.write(header)
        for k, v in line_dict.items():
            fo.write(k + "," + v)
    return outfile


def combine_results(summary_file, samtools_cov_file, pangolin_file, nextclade_file, postrun_output):
    # function that takes in any files generated from this script and adds to the summary dataframe

    print("Generating post-run summary report" + "\n")

    sum_df = pd.read_csv(summary_file, dtype = str)
    sum_df['depth_after_trimming'] = None
    sum_df['1X_coverage_after_trimming'] = None
    columns = list(sum_df.columns)
    columns.remove('Sample')


    if exists(samtools_cov_file):
        print("Getting samtools coverage and depth information from " + samtools_cov_file)
        samtools_df = pd.read_csv(samtools_cov_file, dtype=str, usecols=['Sample', 'coverage', 'meandepth'], index_col=False)
        samtools_df = samtools_df.add_prefix('samtools_')
        st_cols = list(samtools_df.columns)
        st_cols.remove('samtools_Sample')
        st_cols.remove('samtools_meandepth')
        st_cols.remove('samtools_coverage')


        sum_df = pd.merge(sum_df, samtools_df, left_on = "Sample", right_on="samtools_Sample", how = 'outer')
        sum_df['depth_after_trimming'].fillna(sum_df['samtools_meandepth'], inplace=True)
        sum_df['1X_coverage_after_trimming'].fillna(sum_df['samtools_coverage'], inplace=True)
        columns = ['Sample'] + columns
    else:
        print("Combined coverage and depth samtools file wasn't created/doesn't exist.")

    if exists(pangolin_file):
        print("Getting pangolin data from " + pangolin_file)
        p_df = pd.read_csv(pangolin_file, dtype=str)
        p_df = p_df.add_prefix('pangolin_')
        p_cols = list(p_df.columns)
        p_cols.remove('pangolin_taxon')
        p_cols.remove('pangolin_lineage')
    
        sum_df = pd.merge(sum_df, p_df, left_on="Sample", right_on="pangolin_taxon", how = 'outer')
        sum_df.drop('pangolin_taxon', axis=1, inplace=True)
        columns = columns + p_cols
    else:
        print("Combined pangolin file wasn't created/doesn't exist.")

    if exists(nextclade_file):
        print("Getting nextclade data from " + nextclade_file)
        n_df = pd.read_csv(nextclade_file, dtype=str, usecols=['seqName', 'qc.overallScore', 'qc.overallStatus'])
        n_df = n_df.add_prefix('nextclade_')
        n_cols = list(n_df.columns)
        n_cols.remove('nextclade_seqName')
    
        sum_df = pd.merge(sum_df, n_df, left_on="Sample", right_on="nextclade_seqName", how = 'outer')
        columns = columns + n_cols
    else:
        print("Combined nextclade file wasn't created/doesn't exist.")
    
    sum_df.to_csv(postrun_output, columns = columns, index=False)

    # put the combined summary file into an S3 bucket. 
    print("Copying postrun report file to S3 Bucket..." + "\n")
    try:
        with open(postrun_output, "rb") as po:
            s3.put_object(Bucket=outbucket,
            Body=po,
            Key=prefix + "multiqc/" + postrun_output)
        print(f"Successfully copied {postrun_output} to {outbucket}/{prefix}")
    except:
        print(f"Copying {postrun_output} to S3 Bucket {outbucket} was unsuccessful..." + "\n")


def main():
    sum_file = get_summary_file(outdir, "summary_variants_metrics_mqc.csv")
    bam_list = get_bam_files(outdir)
    pangolin_file = get_pangolin_data("int_combined_pangolin.csv", "final_combined_pangolin.csv")
    nextclade_file = get_nextclade_data("int_nextclade.csv", "final_combined_nextclade.csv")
    samtools_file = samtools_coverage("combined_samtools_cov_report.csv", bam_list)
    combine_results(sum_file, samtools_file, pangolin_file, nextclade_file, "postrun_variants_metrics_mqc.csv")


if __name__ == "__main__":
    sys.exit(main())
