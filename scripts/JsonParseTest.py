#!/usr/bin/env python3

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

    # Capture all of the samples index sequence

    indexSequence = []

    for sample in range(len(conversionResults[0]["DemuxResults"])):
        indexSequence.append(conversionResults[0]["DemuxResults"][sample]["IndexMetrics"][0]["IndexSequence"])

    print(indexSequence)

    # Separate out the index 1 and index 2

    index1Sequence = []
    index2Sequence = []

    for index in range(len(indexSequence)):
        splitIndex = indexSequence[index].split("+")
        index1Sequence.append(splitIndex[0])
        index2Sequence.append(splitIndex[1])

    print(index1Sequence)
    print(index2Sequence)

    # Make all possible index1-index2 combinations

    mismatchIndexSequences = []

    for index1 in range(len(index1Sequence)):
        for index2 in range(len(index2Sequence)):
            if index1 != index2:
                mismatchIndexSequences.append(index1Sequence[index1] + "+" + index2Sequence[index2])

    print(mismatchIndexSequences)
    print("Number of mismatched index combinations is " + str(len(mismatchIndexSequences)))

    # Calculate the total number of reads for all samples

    totalNumberOfReads = 0

    for lane in range(numOfLanes):
        for sample in range(numOfSamples):
            numberOfReads = int(conversionResults[lane]["DemuxResults"][sample]["NumberReads"])
            sampleId = conversionResults[lane]["DemuxResults"][sample]["SampleId"]
            #print("For sample " + sampleId + "in lane " + str(lane + 1) + " the number of reads is " + str(numberOfReads))
            totalNumberOfReads += numberOfReads
            #print("The current total number of reads is " + str(totalNumberOfReads))

    # Calculate the mismatched reads

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
                for index1 in index1Sequence:
                    similarity_score1 = compare_sequences(index1, unknown_index1)
                    if similarity_score1 <= 1: # If index 1 is similar or same
                        for index2 in index2Sequence:
                            # Avoid counting a valid index1-index2 combination that is off by 1 base
                            if (index1Sequence.index(index1) != index2Sequence.index(index2)):
                                similarity_score2 = compare_sequences(index2, unknown_index2)
                                if similarity_score2 <= 1: # If index 1 and index 2 are similar or the same
                                    '''
                                    Note: This will not double count the exact same index1-index2 sequence in the
                                    mismatchIndexSequences list because of the if statement above
                                    '''
                                    if barcode not in similar_mismatch_index_dict:
                                        similar_mismatch_index_dict[barcode] = number_of_mismatched_reads
                                    else:
                                        similar_mismatch_index_dict[barcode] += number_of_mismatched_reads
                                else: # If index 1 is same/similar and index 2 isn't
                                    if barcode not in not_similar_mismatch_index_dict:
                                        not_similar_mismatch_index_dict[barcode] = number_of_mismatched_reads
                                    else:
                                        not_similar_mismatch_index_dict[barcode] += number_of_mismatched_reads
                    else: #If index 1 is not similar
                        for index2 in index2Sequence:
                            similarity_score2 = compare_sequences(index2, unknown_index2)
                            if similarity_score2 <= 1: # If index 1 is not similar but index 2 is same/similar
                                if barcode not in not_similar_mismatch_index_dict:
                                    not_similar_mismatch_index_dict[barcode] = number_of_mismatched_reads
                                else:
                                    not_similar_mismatch_index_dict[barcode] += number_of_mismatched_reads
                            else: # If index 1 and index 2 are both not similar, we don't care about this barcode
                                continue

    # Calculate the index hopping percent

    sum_of_mismatched_index_reads = sum(mismatch_index_dict.values())
    sum_of_similar_mismatched_index_reads = sum(similar_mismatch_index_dict.values())

    totalNumberOfMismatchedReads = sum_of_mismatched_index_reads + sum_of_similar_mismatched_index_reads

    indexHopPercent = (totalNumberOfMismatchedReads / totalNumberOfReads) * 100

    print("Similar mismatched Index Dict")
    print(similar_mismatch_index_dict)
    print("Not similar mismatched index dict")
    print(not_similar_mismatch_index_dict)
    print("The total number of mismatched reads is " + str(totalNumberOfMismatchedReads))
    print("The total number of identified reads is " + str(totalNumberOfReads))
    print("Index Hopping Percent is " + str(round(indexHopPercent, 2)) + "%")

    # Calculate how many times each index pair jumped and store it in a dictionary

    indexJumpDict = {}

    for index in range(len(indexSequence)):
        indexJumpCounter = 0
        index1 = index1Sequence[index]
        index2 = index2Sequence[index]
        for keys in mismatch_index_dict.keys():
            if (index1Sequence[index] in keys) or (index2Sequence[index] in keys):
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
        indexJumpDict[indexSequence[index]] = indexJumpCounter

    # Make result file

    resultFile = open(outPath + flowCellId + "_Results.txt", 'w')

    resultFile.write("Number of mismatched reads\t" + str(totalNumberOfMismatchedReads) + "\n")
    resultFile.write("Number of identified reads\t" + str(totalNumberOfReads) + "\n")
    resultFile.write("Index Hopping Percent\t" + str(indexHopPercent) + "%\n\n")
    resultFile.write("Sample\tIndex\tIndex Jump Count\n")

    for (key, val) in indexJumpDict.items():
        sampleListIndex = indexSequence.index(key)
        resultFile.write(sampleList[sampleListIndex] + "\t" + key + "\t" + str(val) + "\n")

    jsonFile.close()
    resultFile.close()

    return indexHopPercent, mismatch_index_dict, indexJumpDict, indexSequence

def import_file(jsonFile):
    openJson = json.load(jsonFile)
    flowCellId = openJson["Flowcell"]
    print(flowCellId)
    conversionResults = openJson["ConversionResults"]
    unknownBarcodes = openJson["UnknownBarcodes"]
    # Calculate the number of lanes on the flowcell
    numOfLanes = len(conversionResults)
    print("Number of lanes on this flowcell is " + str(numOfLanes))
    # List of samples
    sampleList = []
    for sample in range(len(conversionResults[0]["DemuxResults"])):
        sampleList.append(conversionResults[0]["DemuxResults"][sample]["SampleId"])
    print(sampleList)
    # Number of samples on this flowcell
    numOfSamples = len(conversionResults[0]["DemuxResults"])
    print("Number of samples on this flowcell is " + str(numOfSamples))
    return flowCellId, conversionResults, numOfLanes, numOfSamples, unknownBarcodes, sampleList

# Function to compare the similarity between two sequences
def compare_sequences(sequence1, sequence2):
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
