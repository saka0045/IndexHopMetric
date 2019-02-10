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

    # Calculate the total number of mismatched barcodes and store the result in a dictionary

    totalNumberOfMismatchedReads = 0
    misMatchIndexDict = {}

    for lane in range(numOfLanes):
        for index in range(len(mismatchIndexSequences)):
            mismatchIndex = mismatchIndexSequences[index]
            if mismatchIndex in unknownBarcodes[lane]["Barcodes"]:
                numberOfMismatchedReads = int(unknownBarcodes[lane]["Barcodes"][mismatchIndex])
                if mismatchIndex not in misMatchIndexDict:
                    misMatchIndexDict[mismatchIndex] = numberOfMismatchedReads
                else:
                    misMatchIndexDict[mismatchIndex] += numberOfMismatchedReads
                #print("The mismatch index " + mismatchIndex + " was found and has a read number of " + str(numberOfMismatchedReads) + " for lane " + str(lane + 1))
                totalNumberOfMismatchedReads += numberOfMismatchedReads
                #print("The total number of mismatch reads is " + str(totalNumberOfMismatchedReads))
            else:
                misMatchIndexDict[mismatchIndex] = 0
                print("The mismatch index " + mismatchIndex + " was not found")
                #print("The total number of mismatch reads is " + str(totalNumberOfMismatchedReads))

    # Calculate the index hopping percent

    indexHopPercent = (totalNumberOfMismatchedReads / totalNumberOfReads) * 100

    print("The total number of mismatched reads is " + str(totalNumberOfMismatchedReads))
    print("The total number of identified reads is " + str(totalNumberOfReads))
    print("Index Hopping Percent is " + str(round(indexHopPercent, 2)) + "%")

    # Calculate how many times each index pair jumped and store it in a dictionary

    indexJumpDict = {}

    for index in range(len(indexSequence)):
        indexJumpCounter = 0
        for keys in misMatchIndexDict.keys():
            if (index1Sequence[index] in keys) or (index2Sequence[index] in keys):
                #print(index1Sequence[index] + " or " + index2Sequence[index] + " found in " + keys + " with a counter of" + str(misMatchIndexDict[keys]))
                indexJumpCounter += misMatchIndexDict[keys]
        indexJumpDict[indexSequence[index]] = indexJumpCounter
        #print("Total number of index " + indexSequence[index] + " jumped is " + str(indexJumpCounter))

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

    return indexHopPercent, misMatchIndexDict, indexJumpDict, indexSequence

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


if __name__ == "__main__":
    #jsonFile = open("/Users/m006703/IndexHopMetric/files/HH7VKDSXX_Stats.json", 'r')
    #indexHopPercent, misMatchIndexDict, indexJumpDict, indexSequence = main(jsonFile)
    #jsonFile.close()
    main()
