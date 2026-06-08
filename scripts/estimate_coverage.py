#!/bin/env python3

import argparse
import os
import sys
import gzip
from   statistics import mean
from   collections import defaultdict

"""
    Depth columns (0-based)
    0  Chr
    1  BP
    2  DEPTH
"""


def get_arguments() -> tuple:
    """get the arguments"""

    d = "Calculate the coverage statistics from a .depth file"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-d", "--depth", help="Depth file generated from samtools depth", required=True)
    parser.add_argument("-w", "--window_size", help="Size of the window in bp to scan the genome [default = 150,000 bp]", type=int, default=150_000)

    args  = parser.parse_args()
    dfile = args.depth
    wsize = args.window_size

    assert os.path.isfile(dfile), f"Could not locate {dfile}"
    assert wsize >= 0,  "--window_size cannot be less than 0"

    return (dfile, wsize)

def slide_window(depth_file: str, wsize: int) -> None:
    """calculate the average coverage using a sliding window"""

    if (depth_file.endswith(".gz")):
        fh     = gzip.open(depth_file, 'rt')
        oname1 = "GenomeCoverage_" + os.path.basename(depth_file).replace(".depth.gz", '.txt')
        oname2 = oname1.replace("GenomeCoverage_", "WindowDepths_")
        ofh1   = open(oname1, 'w')
        ofh2   = open(oname2, 'w')
    else:
        fh     = open(depth_file, 'r')
        oname1 = "GenomeCoverage_" + os.path.basename(depth_file).replace(".depth", ".txt")
        oname2 = oname1.replace("GenomeCoverage_", "WindowDepths_")
        ofh1   = open(oname1, 'w')
        ofh2   = open(oname2, 'w')

    chr_sum = defaultdict(float)
    chr_cnt = defaultdict(int)
    all_win = defaultdict(int)
    wstart  = 1
    wend    = wsize
    wcnt    = 0
    cur_sum = 0
    cur_chr = ''

    for line in fh:
        fields = line[:-1].split('\t')
        chr    = fields[0]
        bp     = int(fields[1])
        depth  = int(fields[2])
        if (cur_chr == '') and (wstart <= bp <= wend):
            cur_chr = chr
            cur_sum += depth
        elif (chr == cur_chr):
            if (wstart <= bp <= wend):
                cur_sum += depth
            else:
                avg = cur_sum / wsize
                all_win[avg] += 1
                chr_sum[chr] += avg
                chr_cnt[chr] += 1
                wcnt         += 1
                cur_sum = depth
                wstart += wsize
                wend   += wsize
        else:
            avg = cur_sum / wsize
            all_win[avg]     += 1
            chr_sum[cur_chr] += avg
            chr_cnt[cur_chr] += 1
            wcnt             += 1
            cur_sum = depth
            cur_chr = chr
            wstart  = 1
            wend    = wsize

    fh.close()

    # get last bit
    avg = cur_sum / wsize
    all_win[avg]     += 1
    chr_sum[cur_chr] += avg
    chr_cnt[cur_chr] += 1
    wcnt             += 1

    # now to summarize the counts
    for chr in chr_sum.keys():
        sum = chr_sum[chr]
        cnt = chr_cnt[chr]
        cov = round(sum / cnt, 2)
        out = f"{chr} average coverage: {cov}x\n"
        ofh1.write(out)

    ofh1.close()

    header = "WindowAverage\tNumberOfWindows\tPercentOfTotalWindows\n"
    ofh2.write(header)
    for cov, cnt in all_win.items():
        per   = round(cnt/wcnt, 2)
        oline = f"{cov}\t{cnt}\t{per}%\n"
        ofh2.write(oline)
    
    ofh2.close()


def main() -> int:
    """entry point for the small application"""

    # get arguments
    dfile, wsize = get_arguments()

    # parse the depth file and report the coverage
    slide_window(dfile, wsize)

    return 0

if __name__ == "__main__":
    main()
