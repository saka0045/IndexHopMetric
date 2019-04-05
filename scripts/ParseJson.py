#!/usr/bin/env python3
import itertools
import json
import argparse
import sys
import os


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

    args = parser.parse_args()

    inputFile = os.path.abspath(args.inputFile)
    outPath = os.path.abspath(args.outPath)

    # Add / at the end if it is not included in the output path
    if outPath.endswith("/"):
        outPath = outPath
    else:
        outPath = outPath + "/"

    jsonFile = open(inputFile, 'r')

    # Convert the json file to a python dictionary

    flowCellId, conversionResults, numOfLanes, numOfSamples, unknownBarcodes, sampleList = import_file(jsonFile)
    print(flowCellId)
    print(sampleList)
    print("Number of lanes on this flowcell is " + str(numOfLanes))
    print("Number of samples on this flowcell is " + str(numOfSamples))

    # Capture all of the samples index sequence

    indexSequence = capture_index_sequence(conversionResults)
    print(indexSequence)

    # Separate out the index 1 and index 2 and make all possible index1+index2 combinations

    index1Sequence, index2Sequence, mismatchIndexSequences = make_all_possible_index_combinations(indexSequence)
    print(index1Sequence)
    print(index2Sequence)
    print(mismatchIndexSequences)
    print("Number of mismatched index combinations is " + str(len(mismatchIndexSequences)))

    # Calculate the total number of reads for all samples

    totalNumberOfReads = total_number_of_reads(conversionResults, numOfLanes, numOfSamples)

    # Make dictionaries of mismatched reads

    mismatch_index_dict, similar_mismatch_index_dict, not_similar_mismatch_index_dict = mismatched_reads(index1Sequence,
                                                                        index2Sequence, mismatchIndexSequences,
                                                                                    numOfLanes, unknownBarcodes)

    # Calculate the index hopping percent

    indexHopPercent, totalNumberOfMismatchedReads = index_hopping_percent(mismatch_index_dict,
                                                                          similar_mismatch_index_dict,
                                                                          totalNumberOfReads)

    print("Mismatch Index Dict")
    print(mismatch_index_dict)
    print("Similar mismatched Index Dict")
    print(similar_mismatch_index_dict)
    print("Not similar mismatched index dict")
    print(not_similar_mismatch_index_dict)
    print("The total number of mismatched reads is " + str(totalNumberOfMismatchedReads))
    print("The total number of identified reads is " + str(totalNumberOfReads))
    print("Index Hopping Percent is " + str(round(indexHopPercent, 2)) + "%")

    # Calculate how many times each index pair jumped and store it in a dictionary

    indexJumpDict = index_jump_count(index1Sequence, index2Sequence, indexSequence, mismatch_index_dict,
                                     similar_mismatch_index_dict)

    # Group all of the not similar indexes to the valid sample index1+index2 combination

    not_similar_jump_count = {}
    not_similar_index_association = {}
    for index in range(len(indexSequence)):
        not_similar_indexes = {}
        not_similar_index_jump_counter = 0
        index1 = index1Sequence[index]
        index2 = index2Sequence[index]
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
        not_similar_index_association[indexSequence[index]] = not_similar_indexes
        not_similar_jump_count[indexSequence[index]] = not_similar_index_jump_counter

    print("Not similar index association dictionary")
    print(not_similar_index_association)

    # Make result file

    resultFile = open(outPath + flowCellId + "_Results.txt", 'w')

    resultFile.write("Number of mismatched reads\t" + str(totalNumberOfMismatchedReads) + "\n")
    resultFile.write("Number of identified reads\t" + str(totalNumberOfReads) + "\n")
    resultFile.write("Index Hopping Percent\t" + str(indexHopPercent) + "%\n\n")
    resultFile.write("Sample\tIndex\tIndex Jump Count\n")

    for (key, val) in indexJumpDict.items():
        sampleListIndex = indexSequence.index(key)
        resultFile.write(sampleList[sampleListIndex] + "\t" + key + "\t" + str(val) + "\n")

    # Make non similar mismatch index result file

    not_similar_mismatch_index_result_file = open(outPath + flowCellId + "_NonSimilarIndex.txt", 'w')
    not_similar_mismatch_index_result_file.write("Sample\tMismatch Index\tMismatch Read Count\t" +
                                                 "Total Mismatch Reads for Sample\n")

    for (key, val) in not_similar_index_association.items():
        sampleListIndex = indexSequence.index(key)
        not_similar_mismatch_index_result_file.write(sampleList[sampleListIndex] + ": " + key + "\t\t\t"
                                                     + str(not_similar_jump_count[key]) + "\n")
        for (item, number) in val.items():
            not_similar_mismatch_index_result_file.write("\t" + item + "\t" + str(number) + "\n")
        not_similar_mismatch_index_result_file.write("\n")

    jsonFile.close()
    resultFile.close()
    not_similar_mismatch_index_result_file.close()


def index_jump_count(index1Sequence, index2Sequence, indexSequence, mismatch_index_dict, similar_mismatch_index_dict):
    """
    Uses the 3 mismatch index dictionaries to calculate the index jump count, number of how many times a valid
    index1+index2 jumped to a different index1 or index2
    :param index1Sequence:
    :param index2Sequence:
    :param indexSequence:
    :param mismatch_index_dict:
    :param similar_mismatch_index_dict:
    :return Dictionary of the valid index1+index2 combination with the index jump count as values:
    """
    indexJumpDict = {}
    for index_num, index in enumerate(indexSequence):
        indexJumpCounter = 0
        try:
            index1 = index1Sequence[index_num]
        except IndexError:
            raise MissingSequence1Exception
        index2 = index2Sequence[index_num]
        for keys in mismatch_index_dict.keys():
            if (index1Sequence[index_num] in keys) or (index2Sequence[index_num] in keys):
                indexJumpCounter += mismatch_index_dict[keys]
        for keys in similar_mismatch_index_dict.keys():
            sequence1 = keys.split("+")[0]
            sequence2 = keys.split("+")[1]
            similarity_score1 = compare_sequences(index1, sequence1)
            similarity_score2 = compare_sequences(index2, sequence2)
            if similarity_score1 <= 1:
                indexJumpCounter += similar_mismatch_index_dict[keys]
            elif similarity_score2 <= 1:
                indexJumpCounter += similar_mismatch_index_dict[keys]
        '''
        for keys in not_similar_mismatch_index_dict.keys():
            sequence1 = keys.split("+")[0]
            sequence2 = keys.split("+")[1]
            similarity_score1 = compare_sequences(index1, sequence1)
            similarity_score2 = compare_sequences(index2, sequence2)
            if similarity_score1 <= 1:
                indexJumpCounter += not_similar_mismatch_index_dict[keys]
            elif similarity_score2 <= 1:
                indexJumpCounter += not_similar_mismatch_index_dict[keys]
        '''
        indexJumpDict[index] = indexJumpCounter
    return indexJumpDict


def index_hopping_percent(mismatch_index_dict, similar_mismatch_index_dict, totalNumberOfReads):
    """
    Calculates the index hopping percent by dividing the read count totals from mismatch_index_dict and
    similar_mismatch_index_dict by the total number of reads from samples
    :param mismatch_index_dict:
    :param similar_mismatch_index_dict:
    :param totalNumberOfReads:
    :return index hopping percent:
    """
    sum_of_mismatched_index_reads = sum(mismatch_index_dict.values())
    sum_of_similar_mismatched_index_reads = sum(similar_mismatch_index_dict.values())
    totalNumberOfMismatchedReads = sum_of_mismatched_index_reads + sum_of_similar_mismatched_index_reads
    indexHopPercent = (totalNumberOfMismatchedReads / totalNumberOfReads) * 100
    return indexHopPercent, totalNumberOfMismatchedReads


def mismatched_reads(index1Sequence, index2Sequence, mismatchIndexSequences, numOfLanes, unknownBarcodes):
    """
    Takes the index 1 and index 2 sequences to make three dictionaries with read counts
    mismatch_index_dict is a dictionary of the invalid index1+index2 combinations
    similar_mismatch_index_dict is a dictionary of invalid index1+index2 combinations that are off by
    one sequencei that is found in the Barcodes section of the json file
    not_similar_mismatch_index_dict has either the index1 or index2 sequence on one end and a completely unknown
    sequence on the other end in the json file, it will also look for index that are one base off of index1 or index2
    :param index1Sequence:
    :param index2Sequence:
    :param mismatchIndexSequences:
    :param numOfLanes:
    :param unknownBarcodes:
    :return Currently only returns the mismatch_index_dict and similar_mismatch_index_dict
            the not_similar_mismatch_index_dict doesn't have a use at this point:
    """
    mismatch_index_dict = {}
    similar_mismatch_index_dict = {}
    not_similar_mismatch_index_dict = {}
    for lane in range(numOfLanes):
        for barcode in unknownBarcodes[lane]["Barcodes"]:
            number_of_mismatched_reads = int(unknownBarcodes[lane]["Barcodes"][barcode])
            # If the unknown barcode is a combination of mismatched index1-index2, add it straight to the dict
            if barcode in mismatchIndexSequences:
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
                for index1 in index1Sequence:
                    index1_loop += 1
                    similarity_score1 = compare_sequences(index1, unknown_index1)
                    if similarity_score1 <= 1:  # If index 1 is similar or same
                        for index2 in index2Sequence:
                            index2_loop += 1
                            similarity_score2 = compare_sequences(index2, unknown_index2)
                            # Avoid counting a valid index1-index2 combination that is off by 1 base
                            if similarity_score2 <= 1 and index1Sequence.index(index1) == index2Sequence.index(index2):
                                index1_loop = 0
                                break
                            elif (index1Sequence.index(index1) != index2Sequence.index(index2)):
                                if similarity_score2 <= 1:  # If index 1 and index 2 are similar or the same
                                    '''
                                    Note: This will not double count the exact same index1-index2 sequence in the
                                    mismatchIndexSequences list because of the if statement above
                                    '''
                                    if barcode not in similar_mismatch_index_dict:
                                        similar_mismatch_index_dict[barcode] = number_of_mismatched_reads
                                        index1_loop = 0
                                        break
                                    else:
                                        similar_mismatch_index_dict[barcode] += number_of_mismatched_reads
                                        index1_loop = 0
                                        break
                                elif index2_loop == len(index2Sequence):  # If index 1 is same/similar and index 2 isn't
                                    # Loop through all the index 2s first before adding the number with no match
                                    if barcode not in not_similar_mismatch_index_dict:
                                        not_similar_mismatch_index_dict[barcode] = number_of_mismatched_reads
                                        index1_loop = 0
                                    else:
                                        not_similar_mismatch_index_dict[barcode] += number_of_mismatched_reads
                                        index1_loop = 0
                    elif index1_loop == len(index1Sequence):  # If index 1 is not similar
                        # Loop through all the index 1 first before adding the read counts
                        for index2 in index2Sequence:
                            similarity_score2 = compare_sequences(index2, unknown_index2)
                            if similarity_score2 <= 1:  # If index 1 is not similar but index 2 is same/similar
                                if barcode not in not_similar_mismatch_index_dict:
                                    not_similar_mismatch_index_dict[barcode] = number_of_mismatched_reads
                                else:
                                    not_similar_mismatch_index_dict[barcode] += number_of_mismatched_reads
                            else:  # If index 1 and index 2 are both not similar, we don't care about this barcode
                                continue
    return mismatch_index_dict, similar_mismatch_index_dict, not_similar_mismatch_index_dict


def total_number_of_reads(conversionResults, numOfLanes, numOfSamples):
    """
    Calculates the total number of reads assigned to samples
    :param conversionResults:
    :param numOfLanes:
    :param numOfSamples:
    :return total number of reads assigned to samples:
    """
    totalNumberOfReads = 0
    for lane in range(numOfLanes):
        for sample in range(numOfSamples):
            numberOfReads = int(conversionResults[lane]["DemuxResults"][sample]["NumberReads"])
            totalNumberOfReads += numberOfReads
    return totalNumberOfReads


def make_all_possible_index_combinations(indexSequence):
    """
    Separates out the indexes in indexSequence to index1 and index2
    and makes all possible index1+index2 combinations
    Assumes the indexes in indexSequence is dual indexed
    :param indexSequence:
    :return index1Sequence, index2Sequence, mismatchIndexSequences:
    """
    index1Sequence = []
    index2Sequence = []
    for index in range(len(indexSequence)):
        splitIndex = indexSequence[index].split("+")
        index1Sequence.append(splitIndex[0])
        index2Sequence.append(splitIndex[1])
    # Make all possible index1-index2 combinations
    mismatchIndexSequences = []
    for index1 in range(len(index1Sequence)):
        for index2 in range(len(index2Sequence)):
            if index1 != index2:
                mismatchIndexSequences.append(index1Sequence[index1] + "+" + index2Sequence[index2])
    return index1Sequence, index2Sequence, mismatchIndexSequences


def capture_index_sequence(conversionResults):
    """
    Capture the index sequences used on samples in a Stats.json file
    :param conversionResults:
    :return indexSequence:
    """
    indexSequence = []
    for sample in range(len(conversionResults[0]["DemuxResults"])):
        indexSequence.append(conversionResults[0]["DemuxResults"][sample]["IndexMetrics"][0]["IndexSequence"])
    return indexSequence


def import_file(jsonFile):
    """
    Import the Stats.json file and convert it to a dictionary
    :param jsonFile:
    :return flowcell ID, number of lanes used, list of sample names, number of samples,
    path to ConversionResults and UnknownBarcodes:
    """
    openJson = json.load(jsonFile) # move this to the main and added it to a variable to pass along the function
    flowCellId = openJson["Flowcell"]
    conversionResults = openJson["ConversionResults"]
    unknownBarcodes = openJson["UnknownBarcodes"]
    # Calculate the number of lanes on the flowcell
    numOfLanes = len(conversionResults)
    # List of samples
    sampleList = []
    for sample in range(len(conversionResults[0]["DemuxResults"])):
        sampleList.append(conversionResults[0]["DemuxResults"][sample]["SampleId"])
    # Number of samples on this flowcell
    numOfSamples = len(conversionResults[0]["DemuxResults"])
    return flowCellId, conversionResults, numOfLanes, numOfSamples, unknownBarcodes, sampleList


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


if __name__ == "__main__":
    main()


class MissingSequence1Exception(Exception):
    pass
