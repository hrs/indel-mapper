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
    parser.add_argument('-o', '--output', help="Filename to save the output as", required=True)

    return parser.parse_args()

def write_csv(filename, results):
    output_file = open(filename, "w")
    writer = csv.writer(output_file)

    header = ["Name", "Sequence", "N20", "PAM", "Total Reads", "Cutsite", "Count", "Reference", "Read"]
    writer.writerow(header)

    for reference_presenter in results:
        sequence_info = [reference_presenter.name(),
                         reference_presenter.sequence(),
                         reference_presenter.n20(),
                         reference_presenter.pam(),
                         reference_presenter.total_reads()]

        for cluster in reference_presenter.mutation_clusters:
            cluster_info = [cluster.cutsite_region, cluster.count()]
            for sequence in cluster.representations:
                read_info = [sequence.reference, sequence.read]
                writer.writerow(sequence_info + cluster_info + read_info)
            if cluster.count() == 0:
                writer.writerow(sequence_info + cluster_info + ["", ""])
        if len(reference_presenter.mutation_clusters) == 0:
            writer.writerow(sequence_info + ["", "", "", ""])

    output_file.close()

def main(args):

    reference_name_to_reads = SamParser(pysam.AlignmentFile(args.alignment, "rb")).reference_name_to_reads_dict()
    references = ReferenceParser(csv.reader(open(args.reference)), reference_name_to_reads).references()

    presenter_results = Presenter([reference for reference in references if reference.is_valid]).present()

    write_csv(args.output, presenter_results)

args = get_args()

if __name__ == '__main__':
    main(args)
