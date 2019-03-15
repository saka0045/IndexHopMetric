import argparse
import gzip
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-u', '--underminedPath', dest='undetermined_path', required=True,
        help='Path to fastq file'
    )
    parser.add_argument(
        '-o', '--outPath', dest='out_path', required=True,
        help="Output path for result file"
    )

    args = parser.parse_args()

    undetermined_path = os.path.abspath(args.undetermined_path)
    out_path = os.path.abspath(args.out_path)

    # Add / at the end of input and out_path if it is not included

    if undetermined_path.endswith("/"):
        undetermined_path = undetermined_path
    else:
        undetermined_path = undetermined_path + "/"
    if out_path.endswith("/"):
        out_path = out_path
    else:
        out_path = out_path + "/"

    # Make a list of fastq files

    file_list = os.listdir(undetermined_path)

    # Take out the fastq files ad separate out the undetermined fastq files

    fastq_file_list = []
    undetermined_fastq_file_list = []

    for file in file_list:
        if file.endswith("fastq.gz"):
            if file.startswith("Undetermined"):
                undetermined_fastq_file_list.append(file)
            else:
                fastq_file_list.append(file)

    print("Fastq file list")
    print(fastq_file_list)
    print("Undetermined file list")
    print(undetermined_fastq_file_list)

    # Collect barcodes from sample fastq files

    sample_barcode_dict = collect_barcodes(undetermined_path, fastq_file_list)

    # Collect barcodes from undetermined file list

    undetermined_barcode_dict = collect_barcodes(undetermined_path, undetermined_fastq_file_list)

    print(sample_barcode_dict)
    print(undetermined_barcode_dict)


def collect_barcodes(path, filelist):
    barcode_dict = {}

    for file in filelist:
        print("Processing: " + path + file)
        line_number = 0
        fastq_file = gzip.open(path + file, "rt")

        for line in fastq_file:
            if line_number % 4 == 0:
                line = line.strip("\n")
                barcode = line.split(":")[-1]
                if barcode not in barcode_dict:
                    barcode_dict[barcode] = 1
                else:
                    barcode_dict[barcode] += 1
                line_number += 1
            else:
                line_number += 1

    return barcode_dict


if __name__ == "__main__":
    main()
