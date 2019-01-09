#!/usr/bin/env python3

import json
def main(jsonFile):

    # Convert the json file to a python dictionary

    conversionResults, numOfLanes, numOfSamples, unknownBarcodes = import_file(jsonFile)

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
            print("For sample " + sampleId + "in lane " + str(lane + 1) + " the number of reads is " + str(numberOfReads))
            totalNumberOfReads += numberOfReads
            print("The current total number of reads is " + str(totalNumberOfReads))

    # Calculate the total number of mismatched barcodes

    totalNumberOfMismatchedReads = 0

    for lane in range(numOfLanes):
        for mismatchIndex in range(len(mismatchIndexSequences)):
            if mismatchIndexSequences[mismatchIndex] in unknownBarcodes[lane]["Barcodes"]:
                numberOfMismatchedReads = int(unknownBarcodes[lane]["Barcodes"][mismatchIndexSequences[mismatchIndex]])
                print("The mismatch index " + mismatchIndexSequences[mismatchIndex] + " was found and has a read number of " + str(numberOfMismatchedReads) + " for lane " + str(lane + 1))
                totalNumberOfMismatchedReads += numberOfMismatchedReads
                print("The total number of mismatch reads is " + str(totalNumberOfMismatchedReads))
            else:
                print("The mismatch index " + mismatchIndexSequences[mismatchIndex] + " was not found")
                print("The total number of mismatch reads is " + str(totalNumberOfMismatchedReads))

    # Calculate the index hopping percent

    indexHopPercent = (totalNumberOfMismatchedReads / totalNumberOfReads) * 100

    print("The total number of identified reads is " + str(totalNumberOfReads))
    print("Index Hopping Percent is " + str(round(indexHopPercent, 2)) + "%")

    return indexHopPercent


def import_file(jsonFile):
    openJson = json.load(jsonFile)
    conversionResults = openJson["ConversionResults"]
    unknownBarcodes = openJson["UnknownBarcodes"]
    # Calculate the number of lanes on the flowcell
    numOfLanes = len(conversionResults)
    print("Number of lanes on this flowcell is " + str(numOfLanes))
    # Number of samples on this flowcell
    numOfSamples = len(conversionResults[0]["DemuxResults"])
    print("Number of samples on this flowcell is " + str(numOfSamples))
    return conversionResults, numOfLanes, numOfSamples, unknownBarcodes


if __name__ == "__main__":
    jsonFile = open("/Users/m006703/IndexHopMetric/files/Stats.json", 'r')
    main()
