#!/usr/bin/env python3
import argparse
import itertools
import sys

import pysam

def read_one_col(file_path):
    """
    Read one column file and strip the newline character.
    Args:
        file_path: file path.
    Returns:
        list of strings.
    """
    with open(file_path) as fh:
        return [line.strip() for line in fh]


class Barcode:
    """
    ## Features
    - Demultiplex barcodes.
    - Filter invalid R2 reads, which includes:
        - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.
        - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.
        - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
        - Low quality reads: low sequencing quality in barcode and UMI regions.
    ## Output
    - `01.barcode/{sample}_1.fq(.gz)`. Read 1
    - `01.barcode/{sample}_2.fq(.gz)`. Dual index i5 read(Barcode only)
    - `01.barcode/{sample}_3.fq(.gz)`. Read 2
    """

    def __init__(self, args):
        self.args = args
        self.raw_reads = 0
        self.corrected_barcodes = 0
        self.fq1_list = args.fq1.split(',')
        self.fq2_list = args.fq2.split(',')
        self.fq3_list = args.fq3.split(',')

        self.out_fq1 = f'{self.args.sample}_R1.fastq'
        self.out_fq3 = f'{self.args.sample}_R3.fastq'

        self.fh_fq1 = open(self.out_fq1, 'w')
        self.fh_fq3 = open(self.out_fq3, 'w')

    def close_files(self):
        self.fh_fq1.close()
        self.fh_fq3.close()

    @staticmethod
    def get_mismatch_dict(seq_list, n_mismatch=1):
        """
        Return:
        mismatch dict. Key: mismatch seq, value: seq in seq_list

        >>> seq_list = ["AACGTGAT", "AAACATCG"]
        >>> mismatch_dict = Barcode.get_mismatch_dict(seq_list)
        >>> mismatch_dict["AACGTGAA"] == "AACGTGAT"
        True
        """
        mismatch_dict = {}

        for seq in seq_list:
            seq = seq.strip()
            if seq == '':
                continue
            for mismatch_seq in Barcode.findall_mismatch(seq, n_mismatch):
                mismatch_dict[mismatch_seq] = seq

        return mismatch_dict

    @staticmethod
    def parse_whitelist_file(files: list, n_pattern: int, n_mismatch: int):
        """
        files: file paths
        n_pattern: number of sections in pattern
        n_mismatch: allowed number of mismatch bases
        Returns:
            white_set_list
            mismatch_list
        """
        n_files = len(files)
        if n_files == 1 and n_pattern > 1:
            files = [files[0]] * n_pattern
        elif n_files != n_pattern:
            sys.exit(f'number of whitelist files({n_files} files:{files}) != n_pattern({n_pattern})')

        white_set_list, mismatch_list = [], []
        for f in files:
            barcodes, _ = read_one_col(f)
            white_set_list.append(set(barcodes))
            barcode_mismatch_dict = Barcode.get_mismatch_dict(barcodes, n_mismatch)
            mismatch_list.append(barcode_mismatch_dict)

        return white_set_list, mismatch_list

    @staticmethod
    def findall_mismatch(seq, n_mismatch=1, bases='ACGTN'):
        """
        choose locations where there's going to be a mismatch using combinations
        and then construct all satisfying lists using product

        Return:
        all mismatch <= n_mismatch set.

        >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
        >>> seq_set = Barcode.findall_mismatch("ACG")
        >>> seq_set == answer
        True
        """
        seq_set = set()
        seq_len = len(seq)
        if n_mismatch > seq_len:
            n_mismatch = seq_len
        for locs in itertools.combinations(range(seq_len), n_mismatch):
            seq_locs = [[base] for base in seq]
            for loc in locs:
                seq_locs[loc] = list(bases)
            for poss in itertools.product(*seq_locs):
                seq_set.add(''.join(poss))
        return seq_set

    def run(self):
        """
        bool_whitelist = (whitelist_file is not None) and whitelist_file != "None"
        C_len = sum([item[1] - item[0] for item in pattern_dict['C']])

        if bool_whitelist:
            barcode_set_list, barcode_mismatch_list = Barcode.parse_whitelist_file(whitelist_files,
                                            n_pattern=len(pattern_dict['C']), n_mismatch=1)
        """
        for fq1_path,fq2_path,fq3_path in zip(self.fq1_list, self.fq2_list, self.fq3_list):
            with pysam.FastxFile(fq1_path, persist=False) as fq1, \
                    pysam.FastxFile(fq2_path, persist=False) as fq2, \
                        pysam.FastxFile(fq3_path, persist=False) as fq3:


                for entry1, entry2, entry3 in zip(fq1, fq2, fq3):
                    header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    header3, seq3, qual3 = entry3.name, entry3.sequence, entry3.quality
                    self.raw_reads += 1

                    """
                    # barcode filter
                    seq_list = self.get_seq_list(seq2, pattern_dict, 'C')

                    if self.bool_whitelist:
                        bool_valid, bool_corrected, corrected_seq = Barcode.check_seq_mismatch(
                            seq_list, barcode_set_list, barcode_mismatch_list)

                        if not bool_valid:
                            self.no_barcode_num += 1
                            continue
                        elif bool_corrected:
                            self.barcode_corrected_num += 1
                        cb = corrected_seq
                    else:
                        cb = "".join(seq_list)

                    self.clean_num += 1
                    self.barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
                    self.umi_qual_Counter.update(C_U_quals_ascii[C_len:])

                    umi = Barcode.get_seq_str(seq2, pattern_dict['U'])
                    if umi:
                        umi += '_'
                    """
                    cb = seq2
                    umi = ''
                    cb_umi = ':'.join([cb,umi])

                    self.fh_fq1.write(f'@{cb_umi}\n{seq1}\n+\n{qual1}\n')
                    self.fh_fq3.write(f'@{cb_umi}\n{seq3}\n+\n{qual3}\n')

            self.close_files()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', type=str, required=True, help='sample name')
    parser.add_argument('--fq1', type=str, required=True, help='read1 file path')
    parser.add_argument('--fq2', type=str, required=True, help='read2 file path')
    parser.add_argument('--fq3', type=str, required=True, help='read3 file path')
    args = parser.parse_args()

    barcode = Barcode(args)
    barcode.run()
    print(f'raw_reads: {barcode.raw_reads}')
    print(f'corrected_barcodes: {barcode.corrected_barcodes}')
