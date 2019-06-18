#!/usr/bin/env python3

import argparse
import os
from datetime import date

from IndexHoppingLib import parse_sample_sheet, import_file, validate_samples, capture_index_sequence, \
    make_all_possible_index_combinations, calculate_total_number_of_reads, mismatched_reads, index_hopping_percent, \
    index_jump_count, not_similar_jump_count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--inputFile', dest='inputFile', required=True,
        help='Path to input json file'
    )
    parser.add_argument(
        '-o', '--outPath', dest='outPath', required=True,
        help='Output path for result file'
    )
    parser.add_argument(
        '-s', '--sampleSheet', dest='sample_sheet', required=True,
        help="Path to sample sheet"
    )
    parser.add_argument(
        '-r', '--report', action="store_true",
        help="Creates .txt report if argument is present"
    )

    args = parser.parse_args()

    input_file = os.path.abspath(args.inputFile)
    out_path = os.path.abspath(args.outPath)
    sample_sheet_path = os.path.abspath(args.sample_sheet)

    # Add / at the end if it is not included in the output path
    if out_path.endswith("/"):
        out_path = out_path
    else:
        out_path = out_path + "/"

    json_file = open(input_file, 'r')

    # Gather information from Sample Sheet

    sample_sheet_info = parse_sample_sheet(sample_sheet_path)

    # Convert the json file to a python dictionary

    flow_cell_id, run_dir_base, conversion_results, num_of_lanes, num_of_samples, \
    unknown_barcodes, sample_list = import_file(json_file)

    # Validate the sample name from sample list to Stats.json

    validate_samples(sample_list, sample_sheet_info)

    print(flow_cell_id)
    print(sample_list)
    print("Number of lanes on this flowcell is " + str(num_of_lanes))
    print("Number of samples on this flowcell is " + str(num_of_samples))

    # Capture all of the samples index sequence

    index_sequence = capture_index_sequence(conversion_results)
    print(index_sequence)

    # Separate out the index 1 and index 2 and make all possible index1+index2 combinations

    index1_sequence, index2_sequence, mismatch_index_sequences = make_all_possible_index_combinations(index_sequence)
    print(index1_sequence)
    print(index2_sequence)
    # print(mismatch_index_sequences)
    print("Number of mismatched index combinations is " + str(len(mismatch_index_sequences)))

    # Calculate the total number of reads for all samples

    total_number_of_reads = calculate_total_number_of_reads(conversion_results, num_of_lanes, num_of_samples)

    # Make dictionaries of mismatched reads

    mismatch_index_dict, similar_mismatch_index_dict, not_similar_mismatch_index_dict = mismatched_reads(
        index1_sequence,
        index2_sequence, mismatch_index_sequences,
        num_of_lanes, unknown_barcodes)

    # Calculate the index hopping percent

    index_hop_percent, total_number_of_mismatched_reads = index_hopping_percent(mismatch_index_dict,
                                                                                similar_mismatch_index_dict,
                                                                                total_number_of_reads)

    '''
    print("Mismatch Index Dict")
    print(mismatch_index_dict)
    print("Similar mismatched Index Dict")
    print(similar_mismatch_index_dict)
    print("Not similar mismatched index dict")
    print(not_similar_mismatch_index_dict)
    '''
    print("The total number of mismatched reads is " + str(total_number_of_mismatched_reads))
    print("The total number of identified reads is " + str(total_number_of_reads))
    print("Index Hopping Percent is " + str(round(index_hop_percent, 2)) + "%")

    # Calculate how many times each index pair jumped and store it in a dictionary

    index_jump_dict = index_jump_count(index1_sequence, index2_sequence, index_sequence, mismatch_index_dict,
                                       similar_mismatch_index_dict)
    # print(index_jump_dict)

    # Group all of the not similar indexes to the valid sample index1+index2 combination

    not_similar_index_association, not_similar_jump_count_dict = not_similar_jump_count(index1_sequence,
                                                                                        index2_sequence, index_sequence,
                                                                                        not_similar_mismatch_index_dict)

    # print("Not similar index association dictionary")
    # print(not_similar_index_association)

    # Make RoQCM TSV files

    today = date.today().strftime('%Y%m%d')

    index_hopping_tsv = open(out_path + run_dir_base + "-index-hopping-" + today + ".metrictsv", 'w')
    index_hopping_tsv.write("index_hopping\trun.index.hopping.percent\tD\t" + str(index_hop_percent) + "\t"
                            + run_dir_base + "\trun\n")

    index_jump_count_tsv = open(out_path + run_dir_base + "-index-jump-count-" + today + ".metric.tsv", 'w')
    for (key, val) in index_jump_dict.items():
        sample_list_index = index_sequence.index(key)
        sample_id = sample_list[sample_list_index]
        batch_id = sample_sheet_info[sample_id]["Batch_ID"]
        index_jump_count_tsv.write("index_hopping\tsample.index.jump.count\tI\t" + str(val) + "\t" + sample_id +
                                   "\tsample\t" + batch_id + "\tbatch\t" + run_dir_base + "\trun\n")

    index_hopping_tsv.close()
    index_jump_count_tsv.close()

    # Make result file if -r argument is present

    if args.report:
        result_file = open(out_path + flow_cell_id + "_Results.txt", 'w')

        result_file.write("Number of mismatched reads\t" + str(total_number_of_mismatched_reads) + "\n")
        result_file.write("Number of identified reads\t" + str(total_number_of_reads) + "\n")
        result_file.write("Index Hopping Percent\t" + str(index_hop_percent) + "%\n\n")
        result_file.write("Sample\tIndex\tIndex Jump Count\n")

        for (key, val) in index_jump_dict.items():
            sample_list_index = index_sequence.index(key)
            result_file.write(sample_list[sample_list_index] + "\t" + key + "\t" + str(val) + "\n")

        # Make non similar mismatch index result file

        not_similar_mismatch_index_result_file = open(out_path + flow_cell_id + "_NonSimilarIndex.txt", 'w')
        not_similar_mismatch_index_result_file.write("Sample\tMismatch Index\tMismatch Read Count\t" +
                                                     "Total Mismatch Reads for Sample\n")

        for (key, val) in not_similar_index_association.items():
            sample_list_index = index_sequence.index(key)
            not_similar_mismatch_index_result_file.write(sample_list[sample_list_index] + ": " + key + "\t\t\t"
                                                         + str(not_similar_jump_count_dict[key]) + "\n")
            for (item, number) in val.items():
                not_similar_mismatch_index_result_file.write("\t" + item + "\t" + str(number) + "\n")
            not_similar_mismatch_index_result_file.write("\n")

        result_file.close()
        not_similar_mismatch_index_result_file.close()

    json_file.close()


if __name__ == "__main__":
    main()

