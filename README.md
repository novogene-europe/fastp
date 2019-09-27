# fastp

This is a modified version, contains function of length limitation of adapter trimming. For more detailed information, please refer to [fastp](https://github.com/OpenGene/fastp). And a prebuiled binary file based on linux system is under centos directory.



## Release 9/27/2019
1. add --min_trim_length option, which is used to limit the length of adapter trimming, default 10.

2. the origin --overlap_diff_limit option is unapplicable, we revised the code and used it to selection the overlap result. 

3. parse_fastp.py is used to summary the fastp's json file and plot.

For general usage:
```
## only overlap used for trimming
./centos/fastp -i testdata/R1.fq.gz \
  -I testdata/R2.fq.gz \
  -o testdata/R1.clean.fq.gz \
  -O testdata/R2.clean.fq.gz \
  --qualified_quality_phred 5 \
  --unqualified_percent_limit 50 \
  --n_base_limit 15 \
  --overlap_len_require 30 \
  --overlap_diff_limit 1 \
  --min_trim_length 10 \
  --overlap_diff_percent_limit 10 \
  --length_limit 150 \
  --disable_trim_poly_g \
  --json testdata/test_fastp.json \
  --html testdata/test_fastp.html \
  --report_title fastp_test
```

```
##  overlap and adapter sequence trimming
./centos/fastp -i testdata/R1.fq.gz \
  -I testdata/R2.fq.gz \
  -o testdata/R1.clean.fq.gz \
  -O testdata/R2.clean.fq.gz \
  --qualified_quality_phred 5 \
  --unqualified_percent_limit 50 \
  --n_base_limit 15 \
  --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --overlap_len_require 30 \
  --overlap_diff_limit 1 \
  --min_trim_length 10 \
  --overlap_diff_percent_limit 10 \
  --length_limit 150 \
  --disable_trim_poly_g \
  --json testdata/test_fastp.json \
  --html testdata/test_fastp.html \
  --report_title fastp_test
```

## usage: parse_fastp.py

```
Parse fastp output

optional arguments:
  -h, --help            show this help message and exit
  --fastp-json FASTP_JSON, -i FASTP_JSON
                        Fastp json file
  --outdir OUTDIR, -o OUTDIR
                        Output directory
  --samplename SAMPLENAME, -s SAMPLENAME
                        sample name
```
Example:
```
python parse_fastp.py -i testdata/test_fastp.json -o testdata -s test

```
