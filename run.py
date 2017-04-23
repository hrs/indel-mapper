import argparse
import csv
import pysam
from indel_mapper_lib.sam_parser import SamParser
from indel_mapper_lib.reference_parser import ReferenceParser
from indel_mapper_lib.presenter import Presenter

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--alignment', help="Alignment SAM file", required=True)
    parser.add_argument('-r', '--reference', help="Reference CSV file", required=True)
    parser.add_argument('-o', '--output', help="Output file, as CSV", required=True)

    return parser.parse_args()

def write_csv(filename, results):
    output_file = open(filename, "w")
    writer = csv.writer(output_file)
    header = ["Name", "Sequence", "N20", "PAM", "Total Reads", "Cutsite", "Count", "Reference", "Read"]

    writer.writerow(header)
    for reference_presenter in results:
        write_reference(writer, reference_presenter)

    output_file.close()

def write_reference(writer, reference_presenter):
    prefix_cells = reference_presenter.csv_row_prefix_cells()
    if reference_presenter.has_mutation_clusters():
        for cluster in reference_presenter.mutation_clusters:
            cluster_info = [cluster.cutsite_region, cluster.count()]
            for sequence in cluster.representations:
                writer.writerow(prefix_cells + cluster_info + [sequence.reference, sequence.read])
            if cluster.count() == 0:
                writer.writerow(prefix_cells + cluster_info)
    else:
        writer.writerow(prefix_cells)

def main(args):
    reference_name_to_reads = SamParser(pysam.AlignmentFile(args.alignment, "rb")).reference_name_to_reads_dict()
    references = ReferenceParser(csv.reader(open(args.reference)), reference_name_to_reads).references()
    presenter_results = Presenter([reference for reference in references if reference.is_valid]).present()

    write_csv(args.output, presenter_results)

if __name__ == '__main__':
    main(get_args())
