import argparse
import gzip
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--inputPath', dest='input_sample_path', required=True,
        help="Path to the input sample files"
    )
    parser.add_argument(
        '-u', '--undeterminedPath', dest='undetermined_path', required=True,
        help='Path to fastq file'
    )
    parser.add_argument(
        '-o', '--outPath', dest='out_path', required=True,
        help="Output path for result file"
    )
    parser.add_argument(
        '-s', '--sampleSheet', dest='sample_sheet_path', required=True,
        help="Path to sample sheet"
    )

    args = parser.parse_args()

    input_sample_path = os.path.abspath(args.input_sample_path)
    undetermined_path = os.path.abspath(args.undetermined_path)
    out_path = os.path.abspath(args.out_path)
    sample_sheet_path = os.path.abspath(args.sample_sheet_path)

    # Add / at the end of input and out_path if it is not included

    input_sample_path = adjust_path(input_sample_path)
    undetermined_path = adjust_path(undetermined_path)
    out_path = adjust_path(out_path)

    # Make a list of sample fastq files and undetermined fast files

    sample_fastq_file_list = []
    undetermined_fastq_file_list = []

    sample_fastq_file_list, undetermined_fastq_file_list = make_fastq_file_list(input_sample_path,
                                                                                sample_fastq_file_list,
                                                                                undetermined_fastq_file_list)
    sample_fastq_file_list, undetermined_fastq_file_list = make_fastq_file_list(undetermined_path,
                                                                                sample_fastq_file_list,
                                                                                undetermined_fastq_file_list)

    # Print statements for trouble shooting
    print("Fastq file list")
    print(sample_fastq_file_list)
    print("Undetermined file list")
    print(undetermined_fastq_file_list)

    # Collect barcodes from sample fastq files

    sample_barcode_dict = collect_barcodes(input_sample_path, sample_fastq_file_list)

    # Collect barcodes from undetermined file list

    undetermined_barcode_dict = collect_barcodes(undetermined_path, undetermined_fastq_file_list)

    # Calculate reads per sample

    sample_read_count = calculate_sample_read_count(sample_barcode_dict)

    # Calculate read count for undetermined

    undetermined_read_count = calculate_sample_read_count(undetermined_barcode_dict)

    print(sample_barcode_dict)
    for sample in sample_barcode_dict:
        print("Total number of barcodes found in " + sample + ": " + str(sample_read_count[sample]))

    print(undetermined_barcode_dict)
    for sample in undetermined_barcode_dict:
        print("Total number of barcodes found in " + sample + ": " + str(undetermined_read_count[sample]))

    # Parse sample sheet

    print("Opening sample sheet at :" + sample_sheet_path)

    sample_sheet_info = parse_sample_sheet(sample_sheet_path)

    print(sample_sheet_info)


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
                index1 = line_item_list[header_list.index("INDEX")]
                index2 = line_item_list[header_list.index("INDEX2")]
                combined_index = index1 + "+" + index2
                sample_info["Index"] = combined_index
                if sample_id not in sample_sheet_info.keys():
                    sample_sheet_info[sample_id] = sample_info
            break
    sample_sheet.close()

    return sample_sheet_info


def calculate_sample_read_count(sample_barcode_dict):
    """
    Calculates the read count per sample in the respective barcode dictionary with read counts
    :param sample_barcode_dict:
    :return Dictionary of read counts per sample:
    """
    sample_read_count = {}
    for sample in sample_barcode_dict.keys():
        read_count = sum(sample_barcode_dict[sample].values())
        sample_read_count[sample] = read_count

    return sample_read_count


def collect_barcodes(path, file_list):
    """
    Parse the fastq files to collect the barcode information and store it in a dictionary
    The dictionary contains the name of the sample and the associated barcodes it found
    :param path:
    :param file_list:
    :return Dictionary of samples and the associated barcodes:
    """
    barcode_dict = {}

    for file in file_list:
        # Print statement for trouble shooting
        print("Processing: " + path + file)
        sample_name = file.split("_")[0]
        sample_barcode_collection_dict = {}
        line_number = 0
        fastq_file = gzip.open(path + file, "rt")

        for line in fastq_file:
            '''
            First line of every 4 lines in the fastq file contains the metadata
            Assumes the barcode information is the last information separated by ":"
            This format was adopted after Casava 1.8
            Newer version of Illumina software may not have the barcode at the end of the fastq file
            '''
            if line_number % 4 == 0:
                line = line.strip("\n")
                barcode = line.split(":")[-1]
                if barcode not in sample_barcode_collection_dict:
                    sample_barcode_collection_dict[barcode] = 1
                else:
                    sample_barcode_collection_dict[barcode] += 1
                line_number += 1
            else:
                line_number += 1

        if sample_name in barcode_dict:
            for (barcode, read_count) in sample_barcode_collection_dict.items():
                if barcode in barcode_dict[sample_name].keys():
                    barcode_dict[sample_name][barcode] += read_count
                else:
                    barcode_dict[sample_name][barcode] = read_count
        else:
            barcode_dict[sample_name] = sample_barcode_collection_dict

    return barcode_dict


def adjust_path(path):
    """
    Adds a "/" at the end of the directory location if not present for the input paths
    :param path:
    :return returns a path with "/" at the end:
    """
    if path.endswith("/"):
        path = path
    else:
        path = path + "/"

    return path


def make_fastq_file_list(path, sample_fastq_file_list, undetermined_fastq_file_list):
    """
    Makes a list of sample fastqs and undetermined fastqs from the specified path
    :param path:
    :param sample_fastq_file_list:
    :param undetermined_fastq_file_list:
    :return list of sample fastqs and undetermined fastqs:
    """
    file_list = os.listdir(path)

    for file in file_list:
        if file.endswith("fastq.gz"):
            if file.startswith("Undetermined"):
                undetermined_fastq_file_list.append(file)
            else:
                sample_fastq_file_list.append(file)

    return sample_fastq_file_list, undetermined_fastq_file_list


if __name__ == "__main__":
    main()
