#!/usr/bin/env python3
import itertools
import json
import argparse
import sys
import os
from datetime import date


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

    mismatch_index_dict, similar_mismatch_index_dict, not_similar_mismatch_index_dict = mismatched_reads(index1_sequence,
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


def validate_samples(sample_list, sample_sheet_info):
    for sample_id in sample_sheet_info:
        if len(sample_sheet_info) != len(sample_list):
            print("Number of samples between Sample Sheet and Stats.json does not match! "
                  "Check to see if the right Sample Sheet was provided")
            sys.exit()
        if sample_id not in sample_list:
            print(sample_id + " is not found in Stats.json! Check to see if the right Sample Sheet was provided")
            print("Terminating Index Hopping Script")
            sys.exit()


def not_similar_jump_count(index1_sequence, index2_sequence, index_sequence, not_similar_mismatch_index_dict):
    """

    :param index1_sequence:
    :param index2_sequence:
    :param index_sequence:
    :param not_similar_mismatch_index_dict:
    :return: not_similar_index_association dictionary that contains the index1+index2 combinations where only
    1 of the index matched to the valid index1+index2. not_similar_jump_count_dict summarizes the total number
    of invalid index1+index2 counts in the association dictionary to the valid index1+index2 combination
    """
    not_similar_jump_count_dict = {}
    not_similar_index_association = {}
    for index in range(len(index_sequence)):
        not_similar_indexes = {}
        not_similar_index_jump_counter = 0
        index1 = index1_sequence[index]
        index2 = index2_sequence[index]
        for keys in not_similar_mismatch_index_dict.keys():
            sequence1 = keys.split("+")[0]
            sequence2 = keys.split("+")[1]
            similarity_score1 = compare_sequences(index1, sequence1)
            similarity_score2 = compare_sequences(index2, sequence2)
            if similarity_score1 <= 1:
                not_similar_index_jump_counter += not_similar_mismatch_index_dict[keys]
                not_similar_indexes[keys] = not_similar_mismatch_index_dict[keys]
            elif similarity_score2 <= 1:
                not_similar_index_jump_counter += not_similar_mismatch_index_dict[keys]
                not_similar_indexes[keys] = not_similar_mismatch_index_dict[keys]
        not_similar_index_association[index_sequence[index]] = not_similar_indexes
        not_similar_jump_count_dict[index_sequence[index]] = not_similar_index_jump_counter
    return not_similar_index_association, not_similar_jump_count_dict


def index_jump_count(index1_sequence, index2_sequence, index_sequence,
                     mismatch_index_dict, similar_mismatch_index_dict):
    """
    Uses the 3 mismatch index dictionaries to calculate the index jump count, number of how many times a valid
    index1+index2 jumped to a different index1 or index2
    :param index1_sequence:
    :param index2_sequence:
    :param index_sequence:
    :param mismatch_index_dict:
    :param similar_mismatch_index_dict:
    :return Dictionary of the valid index1+index2 combination with the index jump count as values:
    """
    index_jump_dict = {}
    for index_num, index in enumerate(index_sequence):
        index_jump_counter = 0
        try:
            index1 = index1_sequence[index_num]
        except IndexError:
            raise MissingSequence1Exception
        index2 = index2_sequence[index_num]
        for keys in mismatch_index_dict.keys():
            if (index1_sequence[index_num] in keys) or (index2_sequence[index_num] in keys):
                index_jump_counter += mismatch_index_dict[keys]
        for keys in similar_mismatch_index_dict.keys():
            sequence1 = keys.split("+")[0]
            sequence2 = keys.split("+")[1]
            similarity_score1 = compare_sequences(index1, sequence1)
            similarity_score2 = compare_sequences(index2, sequence2)
            if similarity_score1 <= 1:
                index_jump_counter += similar_mismatch_index_dict[keys]
            elif similarity_score2 <= 1:
                index_jump_counter += similar_mismatch_index_dict[keys]
        '''
        for keys in not_similar_mismatch_index_dict.keys():
            sequence1 = keys.split("+")[0]
            sequence2 = keys.split("+")[1]
            similarity_score1 = compare_sequences(index1, sequence1)
            similarity_score2 = compare_sequences(index2, sequence2)
            if similarity_score1 <= 1:
                index_jump_counter += not_similar_mismatch_index_dict[keys]
            elif similarity_score2 <= 1:
                index_jump_counter += not_similar_mismatch_index_dict[keys]
        '''
        index_jump_dict[index] = index_jump_counter
    return index_jump_dict


def index_hopping_percent(mismatch_index_dict, similar_mismatch_index_dict, total_number_of_reads):
    """
    Calculates the index hopping percent by dividing the read count totals from mismatch_index_dict and
    similar_mismatch_index_dict by the total number of reads from samples
    :param mismatch_index_dict:
    :param similar_mismatch_index_dict:
    :param total_number_of_reads:
    :return index hopping percent:
    """
    sum_of_mismatched_index_reads = sum(mismatch_index_dict.values())
    sum_of_similar_mismatched_index_reads = sum(similar_mismatch_index_dict.values())
    total_number_of_mismatched_reads = sum_of_mismatched_index_reads + sum_of_similar_mismatched_index_reads
    index_hop_percent = (total_number_of_mismatched_reads / total_number_of_reads) * 100
    return index_hop_percent, total_number_of_mismatched_reads


def mismatched_reads(index1_sequence, index2_sequence, mismatch_index_sequences, num_of_lanes, unknown_barcodes):
    """
    Takes the index 1 and index 2 sequences to make three dictionaries with read counts
    mismatch_index_dict is a dictionary of the invalid index1+index2 combinations
    similar_mismatch_index_dict is a dictionary of invalid index1+index2 combinations that are off by
    one sequencei that is found in the Barcodes section of the json file
    not_similar_mismatch_index_dict has either the index1 or index2 sequence on one end and a completely unknown
    sequence on the other end in the json file, it will also look for index that are one base off of index1 or index2
    :param index1_sequence:
    :param index2_sequence:
    :param mismatch_index_sequences:
    :param num_of_lanes:
    :param unknown_barcodes:
    :return Currently only returns the mismatch_index_dict and similar_mismatch_index_dict
            the not_similar_mismatch_index_dict doesn't have a use at this point:
    """
    mismatch_index_dict = {}
    similar_mismatch_index_dict = {}
    not_similar_mismatch_index_dict = {}
    for lane in range(num_of_lanes):
        for barcode in unknown_barcodes[lane]["Barcodes"]:
            number_of_mismatched_reads = int(unknown_barcodes[lane]["Barcodes"][barcode])
            # If the unknown barcode is a combination of mismatched index1-index2, add it straight to the dict
            if barcode in mismatch_index_sequences:
                if barcode not in mismatch_index_dict:
                    mismatch_index_dict[barcode] = number_of_mismatched_reads
                else:
                    mismatch_index_dict[barcode] += number_of_mismatched_reads
            # If not a known combinations, see if the index sequence is similar (1 base off) from the known
            # index sequence
            else:
                split_barcode = barcode.split("+")
                unknown_index1 = split_barcode[0]
                unknown_index2 = split_barcode[1]
                index1_loop = 0
                index2_loop = 0
                for index1 in index1_sequence:
                    index1_loop += 1
                    similarity_score1 = compare_sequences(index1, unknown_index1)
                    if similarity_score1 <= 1:  # If index 1 is similar or same
                        for index2 in index2_sequence:
                            index2_loop += 1
                            similarity_score2 = compare_sequences(index2, unknown_index2)
                            # Avoid counting a valid index1-index2 combination that is off by 1 base
                            if similarity_score2 <= 1 and index1_sequence.index(index1) == index2_sequence.index(index2):
                                index1_loop = 0
                                break
                            elif (index1_sequence.index(index1) != index2_sequence.index(index2)):
                                if similarity_score2 <= 1:  # If index 1 and index 2 are similar or the same
                                    '''
                                    Note: This will not double count the exact same index1-index2 sequence in the
                                    mismatch_index_sequences list because of the if statement above
                                    '''
                                    if barcode not in similar_mismatch_index_dict:
                                        similar_mismatch_index_dict[barcode] = number_of_mismatched_reads
                                        index1_loop = 0
                                        break
                                    else:
                                        similar_mismatch_index_dict[barcode] += number_of_mismatched_reads
                                        index1_loop = 0
                                        break
                                elif index2_loop == len(index2_sequence):  # If index 1 is same/similar and index 2 isn't
                                    # Loop through all the index 2s first before adding the number with no match
                                    if barcode not in not_similar_mismatch_index_dict:
                                        not_similar_mismatch_index_dict[barcode] = number_of_mismatched_reads
                                        index1_loop = 0
                                    else:
                                        not_similar_mismatch_index_dict[barcode] += number_of_mismatched_reads
                                        index1_loop = 0
                    elif index1_loop == len(index1_sequence):  # If index 1 is not similar
                        # Loop through all the index 1 first before adding the read counts
                        for index2 in index2_sequence:
                            similarity_score2 = compare_sequences(index2, unknown_index2)
                            if similarity_score2 <= 1:  # If index 1 is not similar but index 2 is same/similar
                                if barcode not in not_similar_mismatch_index_dict:
                                    not_similar_mismatch_index_dict[barcode] = number_of_mismatched_reads
                                else:
                                    not_similar_mismatch_index_dict[barcode] += number_of_mismatched_reads
                            else:  # If index 1 and index 2 are both not similar, we don't care about this barcode
                                continue
    return mismatch_index_dict, similar_mismatch_index_dict, not_similar_mismatch_index_dict


def calculate_total_number_of_reads(conversionResults, numOfLanes, numOfSamples):
    """
    Calculates the total number of reads assigned to samples
    :param conversionResults:
    :param numOfLanes:
    :param numOfSamples:
    :return total number of reads assigned to samples:
    """
    total_number_of_reads = 0
    for lane in range(numOfLanes):
        for sample in range(numOfSamples):
            numberOfReads = int(conversionResults[lane]["DemuxResults"][sample]["NumberReads"])
            total_number_of_reads += numberOfReads
    return total_number_of_reads


def make_all_possible_index_combinations(indexSequence):
    """
    Separates out the indexes in indexSequence to index1 and index2
    and makes all possible index1+index2 combinations
    Assumes the indexes in indexSequence is dual indexed
    :param indexSequence:
    :return index1_sequence, index2_sequence, mismatch_index_sequences:
    """
    index1_sequence = []
    index2_sequence = []
    for index in range(len(indexSequence)):
        split_index = indexSequence[index].split("+")
        index1_sequence.append(split_index[0])
        index2_sequence.append(split_index[1])
    # Make all possible index1-index2 combinations
    mismatch_index_sequences = []
    for index1 in range(len(index1_sequence)):
        for index2 in range(len(index2_sequence)):
            if index1 != index2:
                mismatch_index_sequences.append(index1_sequence[index1] + "+" + index2_sequence[index2])
    return index1_sequence, index2_sequence, mismatch_index_sequences


def capture_index_sequence(conversion_results):
    """
    Capture the index sequences used on samples in a Stats.json file
    :param conversion_results:
    :return index_sequence:
    """
    index_sequence = []
    for sample in range(len(conversion_results[0]["DemuxResults"])):
        index_sequence.append(conversion_results[0]["DemuxResults"][sample]["IndexMetrics"][0]["IndexSequence"])
    return index_sequence


def import_file(jsonFile):
    """
    Import the Stats.json file and convert it to a dictionary
    :param jsonFile:
    :return flowcell ID, number of lanes used, list of sample names, number of samples,
    path to ConversionResults and UnknownBarcodes:
    """
    open_json = json.load(jsonFile) # move this to the main and added it to a variable to pass along the function
    flow_cell_id = open_json["Flowcell"]
    run_dir_base = open_json["RunId"]
    conversion_results = open_json["ConversionResults"]
    unknown_barcodes = open_json["UnknownBarcodes"]
    # Calculate the number of lanes on the flowcell
    num_of_lanes = len(conversion_results)
    # List of samples
    sample_list = []
    for sample in range(len(conversion_results[0]["DemuxResults"])):
        sample_list.append(conversion_results[0]["DemuxResults"][sample]["SampleId"])
    # Number of samples on this flowcell
    num_of_samples = len(conversion_results[0]["DemuxResults"])
    return flow_cell_id, run_dir_base, conversion_results, num_of_lanes, num_of_samples, unknown_barcodes, sample_list


def compare_sequences(sequence1, sequence2):
    """
    Compares the similarity between sequence 1 and sequence 2
    The two sequences must be of the same length
    :param sequence1:
    :param sequence2:
    :return similarity_score: If the two sequences are identical, the similarity score would be 0
    """
    if (len(sequence1) == len(sequence2)):
        similarity_score = 0
        for nucleotide in range(len(sequence1)):
            if sequence1[nucleotide] != sequence2[nucleotide]:
                similarity_score += 1
        return similarity_score
    else:
        print(sequence1 + " and " + sequence2 + " must be the same length.")
        return


def parse_sample_sheet(sample_sheet_path):
    """
    Opens the supplied sample sheet and parses out information into a dictionary
    :param sample_sheet_path:
    :return Dictionary containing information from sample sheet:
    """
    sample_sheet = open(sample_sheet_path, 'r')
    sample_sheet_info = {}
    for line in sample_sheet:
        # Go to the line with [Data]
        if line.startswith("[Data]"):
            # Skip to the header line
            line = sample_sheet.readline()
            line = line.rstrip()
            header_list = line.split(",")
            for line in sample_sheet:
                sample_info = {}
                line = line.rstrip()
                line_item_list = line.split(",")
                sample_id = line_item_list[header_list.index("SAMPLE_ID")]
                sample_project = line_item_list[header_list.index("SAMPLE_PROJECT")]
                batch_id = sample_project.split("_")[1]
                sample_info["Sample_Project"] = sample_project
                sample_info["Batch_ID"] = batch_id
                if sample_id not in sample_sheet_info.keys():
                    sample_sheet_info[sample_id] = sample_info
            break
    sample_sheet.close()
    return sample_sheet_info


if __name__ == "__main__":
    main()


class MissingSequence1Exception(Exception):
    pass
