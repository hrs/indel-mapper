import argparse
import csv
import pysam
from indel_mapper_lib.sam_parser import SamParser
from indel_mapper_lib.reference_parser import ReferenceParser
from indel_mapper_lib.presenter import Presenter
from indel_mapper_lib.tabular_reference_exporter import TabularReferenceExporter

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--alignment', help="Alignment SAM file", required=True)
    parser.add_argument('-r', '--reference', help="Reference CSV file", required=True)
    parser.add_argument('-o', '--output', help="Output file, as CSV", required=True)

    return parser.parse_args()

def write_csv(filename, results):
    exported_rows = TabularReferenceExporter(results).export()

    with open(filename, "w") as output_file:
        writer = csv.writer(output_file)
        for row in exported_rows:
            writer.writerow(row)

def main(args):
    reference_name_to_reads = SamParser(pysam.AlignmentFile(args.alignment, "rb")).reference_name_to_reads_dict()
    references = ReferenceParser(csv.reader(open(args.reference)), reference_name_to_reads).references()
    presenter_results = Presenter([reference for reference in references if reference.is_valid]).present()

    write_csv(args.output, presenter_results)

if __name__ == '__main__':
    main(get_args())
