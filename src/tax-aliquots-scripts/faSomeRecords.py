#!/usr/bin/env python3
# faSomeRecords.py

import argparse, textwrap
import sys
import gzip
from re import match
from random import choice
from collections import OrderedDict

def write_to_file(store, outfile, found, requested, heads):
    '''
    Function to write fasta to file and print summary message
    to screen.
    '''
    if len(store) > 0:
        with open(outfile, "w") as o:
            for head in heads:
                o.write(">"+head+"\n"+store[head])
        if found == 0:
            print("No sequences found")
        else:
            print(f"Found {found} sequence(s)")
            if found > requested:
                print(f"Found {found-requested} sequence(s) more than requested")
            elif requested > found:
                if requested-found > 0:
                    print(f"Could not find {requested-found} sequence(s)")
                    #print("\n".join(not_found))
            print(f"Sequences saved to: {outfile}")
    else:
        if found == 0:
            print("No sequences found")
        else:
            print(f"Found {found} sequence(s)")
            if found > requested:
                print(f"Found {found-requested} sequence(s) more than requested")
            elif requested > found:
                if requested-found > 0:
                    print(f"Could not find {requested-found} sequence(s)")
                    #print("\n".join(not_found))
            print(f"Sequences saved to: {outfile}")

def parse_fasta(handle, stdout, outfile, order, exclude, requested, joinheads, heads):
    '''
    Main FASTA parser
    '''
    # start seq counter
    found = 0
    # count each record
    total = 0
    # not found list
    not_found = []
    # start ordered dictionary
    store = OrderedDict()
    # check if list order is requested
    if order:
        if exclude:
            for line in handle:
                # checks for header
                if line[0] == ">":
                    # count all
                    total += 1
                    if "|" + line[1:].rstrip().split(" ")[0] + "|" in joinheads:
                        seq = 0
                    else:
                        h = line[1:].rstrip().split(" ")[0]
                        seq = 1
                        store[h] = ''
                        found += 1
                else:
                    if seq == 1:
                        store[h] += line
            requested = total - requested
        else:
            for line in handle:
                # checks for header
                if line[0] == ">":
                    # count all
                    total += 1
                    if "|" + line[1:].rstrip().split(" ")[0] + "|" in joinheads:
                        h = line[1:].rstrip().split(" ")[0]
                        seq = 1
                        store[h] = ''
                        found += 1
                    else:
                        seq = 0
                else:
                    if seq == 1:
                        store[h] += line
        # go through header list and print
        if stdout:
            for head in heads:
                print(">"+head+"\n"+store[head], end='')
        else:
            # write to file
            write_to_file(store, outfile, found, requested, heads)
    else:
        # check if output is to stdout
        if stdout:
            # if the function is to exclude rather than include
            if exclude:
                # it basically inverts the function if exclude is enabled
                for line in handle:
                    # checks for header
                    if line[0] == ">":
                        if "|" + line[1:].rstrip().split(" ")[0] + "|" in joinheads:
                            # just count matches if found
                            seq = 0
                        else:
                            # print everything else
                            seq = 1
                            print(line, end='')
                            found += 1
                    else:
                        # writes the sequence portion
                        if seq == 1:
                            print(line, end='')
            else:
                for line in handle:
                    # checks for header
                    if line[0] == ">":
                        if "|" + line[1:].rstrip().split(" ")[0] + "|" in joinheads:
                            # start writing to screen if found
                            seq = 1
                            print(line, end='')
                            found += 1
                        else:
                            # count missing if no match
                            # not_found += [line[1:].rstrip()]
                            seq = 0
                    else:
                        if seq == 1:
                            print(line, end='')
        else:
            with open(outfile, "w") as o:
                if exclude:
                    # it basically inverts the function if exclude is enabled
                    for line in handle:
                        # checks for header
                        if line[0] == ">":
                            # count all
                            total += 1
                            if "|" + line[1:].rstrip() + "|" in joinheads:
                                # just count matches if found
                                seq = 0
                            else:
                                # print everything else
                                seq = 1
                                o.write(line)
                                found += 1
                        else:
                            # writes the sequence portion
                            if seq == 1:
                                o.write(line)
                    requested = total - requested
                else:
                    for line in handle:
                        # checks for header
                        if line[0] == ">":
                            total += 1
                            if "|" + line[1:].rstrip() + "|" in joinheads:
                                # start writing to screen if found
                                seq = 1
                                o.write(line)
                                found += 1
                            else:
                                # count missing if no match
                                # not_found += [line[1:].rstrip()]
                                seq = 0
                        else:
                            if seq == 1:
                                o.write(line)
                write_to_file(store, outfile, found, requested, heads)

# main function
def main(fasta,records,stdout=False):
    # parse arguments
    parser = argparse.ArgumentParser(prog="faSomeRecords.py",
        formatter_class=argparse.RawTextHelpFormatter,
        description="Retrieves some FASTA records provided a FASTA file and a list or records.",
        epilog=textwrap.dedent('''\
    Records can be sepcified providing a list file with --list/-l
    or one by one in the command line using --records/-r.

    All sequences are now kept in the same input order. So the --keep flag no longer applies

    The elements in the list must be an exact match as those in the
    master FASTA file, and may or may not include the ">" character at
    the begining of the header.

    Example of list file:
    sequence1
    sequence2
    sequence3
    ...
    sequenceN
    '''))

    # check input files and arguments
    # get all headers from input
    heads = records
    if heads[0][0] == ">":
        heads = [ h[1:].split(" ")[0] for h in heads ]

    # get number of requested records
    requested = len(heads)
    # prepare regex for header lookup
    joinheads = "|" + "|".join(heads) + "|"

    # main function
    with open(fasta) as f:
        parse_fasta(f, stdout, "outfile", True, False, requested, joinheads, heads)

if __name__ == "__main__":
    main()
