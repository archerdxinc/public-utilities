#!/usr/bin/env python
""" Subsets a fastq file based on random number or fraction of total"""
from __future__ import division

__author__ = "Aaron Berlin"
__copyright__ = "Copyright 2015, ArcherDX"
__credits__ = ["Aaron Berlin"]
__version__ = "0.1"
__maintainer__ = "Aaron Berlin"
__email__ = "aberlin@archerdx.com"
__status__ = "Production"
__mod_name__ = "subset_fastq_file.py"

import argparse
import sys
import re
import os
from random import random as rand
from itertools import izip
from functools import partial

try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write("Cannot find the Biopython package.  Please install with \"pip install biopython\"\n")
    sys.exit()

DUMMY_PHRED_SCORE = 33

def parse_cmdline_params(arg_list=None):
    """Parses commandline arguments.
    :param arg_list: arguments from the command line
    :type arg_list: list
    :return: dictionary of options
    """

    description = "Create subset of sequence files"

    # Create instance of ArgumentParser
    parser = argparse.ArgumentParser(description=description)

    selection = parser.add_mutually_exclusive_group()

    parser.add_argument("-i",
                        "--input_file",
                        type=argparse.FileType('r'),
                        help="Input Sequence file",
                        default=argparse.SUPPRESS,
                        required=True)

    parser.add_argument("-2",
                        "--read2_input_file",
                        type=argparse.FileType('r'),
                        help="Read 2 Input Sequence file",
                        default=None,
                        required=False)

    parser.add_argument("-o",
                        "--output_header",
                        type=str,
                        help="Output Sequence file header; used as outfile name, with extension appended",
                        required=True)

    selection.add_argument("-R",
                           "--random",
                           type=int,
                           help="Number of sequences to RANDOMLY select",
                           default=None,
                           required=False)

    selection.add_argument("-F",
                           "--fraction",
                           type=float,
                           help="Fraction of total sequences to RANDOMLY select",
                           default=None,
                           required=False)

    # Parse options
    opts = parser.parse_args(args=arg_list)

    opts.include = None
    opts.exclude = None

    if opts.input_file:
        opts.input_file.close()
        opts.input_file = opts.input_file.name

    if opts.read2_input_file:
        opts.read2_input_file.close()
        opts.read2_input_file = opts.read2_input_file.name
    return opts


def log_message(module, message, status=None, file_handle=sys.stderr):
    """
    :param str module: Module calling the logging
    :param str message: Message to log
    :param str status: Status to log - optional
    :param file file_handle: File to log too - optional
    """
    file_handle.write("{} [{}] > {}\n".format(status, module, message))
    file_handle.flush()

log_info = partial(log_message, status='INFO')
log_warning = partial(log_message, status='WARNING')
log_status = partial(log_message, status='STATUS')
log_error = partial(log_message, status='ERROR')


def get_file_type(file_name, strict=True):
    """ Returns the expected file format based on the extension
    :param str file_name: name of file to process
    :param bool strict: Only allow fasta and fastq file types
    :return: File format
    :rtype: str
    """

    # get the file type according to the extension
    # if no extension indicator is present in the filename, default to "txt"
    file_type = os.path.splitext(file_name)[1]
    if len(file_type) == 0:
        return "txt"
    else:
        # strip the leading period from os.path.splitext()'s returned extension
        file_type = file_type[1:]

    # handle file type name aliases/abbreviations
    if file_type == "fa":
        file_type = "fasta"
    if file_type == "fq":
        file_type = "fastq"

    # log an error in the case in which the "strict" flag is set & the file type is neither "fasta" nor "fastq"
    if not (file_type == "fastq" or file_type == "fasta") and strict:
        log_error(__mod_name__, "{} is not a valid file type for this script".format(file_type))
        sys.exit(-1)

    # the return should always be either "fasta", "fastq", or "txt"
    # if the "strict" flag is set, return should be only "fasta" or "fastq"
    return file_type


def get_list_from_file(selection_file):
    """ Gets a list of reads from a text or sequence file
    :param iterable selection_file: file to process
    :return: Selected read names
    :rtype: set
    """
    # set of sequence/record/read? keys to build up
    selection_set = set()
    # allow "fasta", "fastq", or "txt" file types
    file_type = get_file_type(selection_file.name, strict=False)
    # if the input is a fasta or fastq file, the unique reads are unique keys in the dict representation of the file
    if re.search("fast[a|q]", file_type):
        seq_list = SeqIO.to_dict(SeqIO.parse(selection_file.name, file_type))
        selection_set = set(seq_list.keys())
    # if the input is just a basic txt file, the selection is the set of unique lines in the file
    else:
        for line in selection_file:
            selection_set.add(line.rstrip())
    # return a set of the unique reads from the input file
    return selection_set


def select_from_list(in_seqs, selection_list, include=True):
    """ Select reads based on a list
    :param iterable in_seqs: Names to select from
    :param set selection_list: Names to select
    :param bool include: Either include or exclude the names in selection_file
    :return: Selected reads
    :rtype: set
    """
    selected_seqs = set()

    for seq in in_seqs:
        # Inclusion list
        if seq in selection_list and include:
            selected_seqs.add(seq)
            # Found them all - stop looking
            if len(selection_list) == len(selected_seqs):
                break

        # Exclusion list
        if seq not in selection_list and not include:
            selected_seqs.add(seq)

    return selected_seqs


def select_at_random(in_seqs, target_read_count, total_reads=None, logging_fh=sys.stderr):
    """ Randomly selects N reads from a list
    :param iterable in_seqs: Read names to select from
    :param int target_read_count: Number of reads to return
    :param int total_reads: Total reads to select from
    :param str logging_fh: File handle for logging -- optional
    :return: Randomly selected reads
    :rtype: set
    """

    # handles case of float-valued desired read count
    target_read_count = int(target_read_count)

    # This will reduce the efficiency of using a generator by converting it to a set.
    if total_reads is None:
        total_reads = len(list(in_seqs))

    # set of randomly selected sequences to build up
    selected_seqs = set()

    # Handle if we want more reads than provided
    if total_reads < target_read_count:
        log_warning(__mod_name__, "Requested more reads than the file contains, {} reads short of target. "
                                  "Writing all reads" .format(target_read_count - total_reads),
                    file_handle=logging_fh)
        # perform no subsetting if more reads are desired than are available
        return set(in_seqs)

    # If we want exactly the number of reads as in the input just return the full set
    if total_reads == target_read_count:
        return set(in_seqs)

    # running count of the number of reads selected
    reads_found = 0
    # in the expected case (#(desired) < #(available)), iterate over the input sequences available from which to select
    for seq_num, seq_id in enumerate(in_seqs):

        # Check if the target read count has been attained
        if reads_found == target_read_count:
            break

        if rand() < ((target_read_count - reads_found) / (total_reads - seq_num)):

            selected_seqs.add(seq_id)
            reads_found += 1

    return selected_seqs


def generate_seq_list(input_file, file_type):
    """ Load Sequence file into a generator
    :param str input_file: input file
    :param str file_type: Sequencing file type
    :return: Sequence Names
    :rtype: generator
    """
    if input_file is None:
        yield None
    else:
        for record in SeqIO.parse(input_file, file_type):
            yield record.id


def write_sequences_to_file(selected_reads, r1_seqs, output_header, out_type="fastq", r2_seqs=None):
    """ Write sequences to a file
    :param set selected_reads: Sequences selected to be written
    :param iterable r1_seqs: Raw read 1 sequences
    :param str output_header: output file name
    :param str out_type: output file format (fasta, fastq) - optional
    :param iterable r2_seqs: Raw sequences for read 2 - optional
    """
    r2_output = None

    # Prevent double suffix --> safeguard against outfile name that includes file extension
    if output_header.endswith(out_type):
        output_header = re.sub("." + out_type, "", output_header)

    # construct Read 1's output file's fully specified name
    r1_output_file = "{}.{}".format(output_header, out_type)

    if r2_seqs:
        # check is a certain string is in the Read 1 output file name
        if re.search("_R1[.|_]", output_header):
            # if so, fully specify Read 2's output file by swapping R1 with R2
            r2_output_file = "{}.{}".format(re.sub("_R1", "_R2", output_header), out_type)
        # if an R1 indicator is not present in the Read 1 output file name,
        # glue a read number string indicator between the output file names and extension
        else:
            r1_output_file = "{}_R1.{}".format(output_header, out_type)
            r2_output_file = "{}_R2.{}".format(output_header, out_type)

        r2_output = open(r2_output_file, 'w')

    with open(r1_output_file, 'w') as r1_output:
        if r2_seqs:
            for read1, read2 in izip(r1_seqs, r2_seqs):
                if read1.id in selected_reads:
                    SeqIO.write(read1, r1_output, out_type)
                    SeqIO.write(read2, r2_output, out_type)

        else:
            for read1 in r1_seqs:
                if read1.id in selected_reads:
                    SeqIO.write(read1, r1_output, out_type)

    # r2 output file is only opened if r2 sequences are specified
    if r2_seqs:
        r2_output.close()


def generate_output_sequences(input_file, file_type, include=None, exclude=None, random=None, fraction=None,
                              logging_fh=sys.stderr, metrics_file=None):
    """
    Generates the sequence subset
    :param str input_file: Fastq file to subset
    :param str file_type: Type of input_file
    :param list include: Reads to include -- optional
    :param list exclude: Reads to exclude -- optional
    :param int random: Number of reads to include -- optional
    :param float fraction: Fraction of total reads to include -- optional
    :param str logging_fh: File handle for logging -- optional
    :param str metrics_file: File name for metrics -- optional
    :return: Sequences to write out
    :rtype: set
    """
    output_seqs = []

    # count up all the reads in the file
    total_reads = sum(1 for _ in generate_seq_list(input_file, file_type))

    if include:
        log_info(__mod_name__, "Selecting reads for inclusion", file_handle=logging_fh)
        output_seqs = select_from_list(generate_seq_list(input_file, file_type), get_list_from_file(include))
    elif exclude:
        log_info(__mod_name__, "Selecting reads for exclusion", file_handle=logging_fh)
        # override the default flag for inclusion
        output_seqs = select_from_list(generate_seq_list(input_file, file_type), get_list_from_file(exclude), False)
    elif random:
        log_info(__mod_name__, "Selecting {r} reads from {t} at random".format(r=random, t=total_reads),
                 file_handle=logging_fh)
        # random is a number indicating the number of sequences to randomly select
        output_seqs = select_at_random(generate_seq_list(input_file, file_type), random, total_reads)
    elif fraction:
        # handle the passing of a nonsense value for proportion of reads desired
        if 0 < fraction <= 1:
            fraction_as_reads = total_reads * fraction
        else:
            log_error(__mod_name__, "Fraction provided must be between 0 and 1", file_handle=logging_fh)
            sys.exit(-1)
        log_info(__mod_name__, "Selecting {f}% of {t} reads at random".format(f=round(100 * fraction, 4),
                                                                              t=total_reads),
                 file_handle=logging_fh)
        # select the int-truncated (floored) number of reads indicated by the fraction given the total
        output_seqs = select_at_random(generate_seq_list(input_file, file_type), int(fraction_as_reads), total_reads)
    else:
        log_error(__mod_name__, "-R or -F must be specified", file_handle=logging_fh)
        sys.exit()

    # Write out metrics for master workflow
    if metrics_file:
        with open(metrics_file, 'w') as metrics:
            metrics.write("SEQUENCED_READ_COUNT\t{}\n".format(total_reads))
            metrics.write("NORMALIZED_READ_COUNT\t{}\n".format(len(output_seqs)))

    return output_seqs


def workflow(input_file, output_header, read2_input_file=None, include=None, exclude=None, random=None, fraction=None,
             output_format='fastq', logging_fh=sys.stderr, metrics_file=None):
    """
    Main workflow
    :param str input_file: Fastq file to subset
    :param str output_header: Output file name
    :param str read2_input_file: Read 2 fastq file to subset -- optional
    :param list include: Reads to include -- optional
    :param list exclude: Reads to exclude -- optional
    :param int random: Number of reads to include -- optional
    :param float fraction: Fraction of total reads to include -- optional
    :param str output_format: Format for output files -- optional
    :param str logging_fh: File handle for logging -- optional
    :param str metrics_file: File name for metrics -- optional
    """

    # this should always be either "fasta", "fastq", or "txt"
    # just "fasta" or "fastq" if "strict" flag is set in get_file_type() (default)
    file_type = get_file_type(input_file)

    # Log & fail on impossible fasta --> fastq input --> output attempt
    if file_type == 'fasta' and output_format == 'fastq':
        log_warning(__mod_name__, "Fasta input files can't be output as fastq files")
        sys.exit(-1)

    # attempt to generate a subset of the input sequences
    output_seqs = generate_output_sequences(input_file, file_type, include, exclude, random, fraction, logging_fh,
                                            metrics_file)

    # check that sequences were in fact matched; if so, write them to a file
    if output_seqs:
        r2_seqs = None
        
        # If reads are paired end, grab the Reverse sequences
        if read2_input_file is not None:
            r2_seqs = SeqIO.parse(read2_input_file, file_type)
             
        write_sequences_to_file(output_seqs,
                                SeqIO.parse(input_file, file_type),
                                output_header,
                                out_type=output_format,
                                r2_seqs=r2_seqs)
    else:
        log_error(__mod_name__, "No reads were found", file_handle=logging_fh)


def main(args):  # pragma: no cover
    """  Main logic for script
    :param list args: arguments from command line
    :return: prints to stdout
    """

    log_status(__mod_name__, 'Starting')

    # Process command line args -- returns dictionary mapping cmdline arg names to their respective values
    opts = parse_cmdline_params(args[1:])
    workflow(opts.input_file, opts.output_header, opts.read2_input_file, opts.include,
             opts.exclude, opts.random, opts.fraction)

    log_status(__mod_name__, 'Completed')


if __name__ == "__main__":
    main(sys.argv)
