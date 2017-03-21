import pysam

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sam_parser import SamParser
from reference import Reference

def parse_sam(filename):
    reads = SamParser(pysam.AlignmentFile(filename, "rb")).reads()

    reads_with_no_reference = {}
    reference_to_reads = {}

    for read in reads:
        try:
            if read.reference_name in reference_to_reads and len(read.get_reference_positions(full_length=True)) != 0:
                reference_to_reads[read.reference_name].append(read)
            else:
                reference_to_reads[read.reference_name] = [read]
        except:
            reads_with_no_reference[read.query_name] = True

    return reference_to_reads

def get_indels_from_reads(reads, template):

    # has_other means it has another mismatch that's not insertion or deletion
    has_other = False
    reads_with_other_mismatch = []

    all_found_indels = []

    reads_to_indels = {}

    for read in reads:
        reference_positions = read.get_reference_positions(full_length=True)
        indels = get_indels_from_reference_positions(reference_positions)
        has_other = has_other_mismatch(read.cigartuples)

        if has_other:
            has_other_for_reference = True
            reads_with_other_mismatch.append(read.query_name)
        else:
            set_distance_to_cutsite(indels, template)
            all_found_indels.extend(indels)
            reads_to_indels[read.query_name] = indels
    return all_found_indels

def has_other_mismatch(cigartuples):
    for cig in cigartuples:
        cig_type, _ = cig
        if is_other(cig_type):
            return True
    return False

def get_indels_from_reference_positions(reference_positions):

    indels = []

    start_index = reference_positions[0]
    length = 0
    prev_value = reference_positions[0]

    for idx in reference_positions[1:]:
        if idx is not None:
            if prev_value is not None:
                # deletion
                if idx != prev_value + 1 and idx != prev_value:
                    indels.append(Indel(prev_value, idx, (idx-1) - prev_value, True))
            else:
                # end of insertion
                indels.append(Indel(start_index, idx, length, False))
                length = 0
        else:
            # start of insertion
            if prev_value is not None:
                start_index = prev_value
                length += 1
            else:
                length += 1

        prev_value = idx

    return indels

def is_insertion_or_deletion(cigar_type):
    return cigar_type == 1 or cigar_type == 2

def is_other(cigar_type):
    return cigar_type != 0 and not is_insertion_or_deletion(cigar_type)

def set_distance_to_cutsite(indels, template):
    for i in indels:
        i.distance_to_cutsite = template.distance_to_cutsite(i)

def get_reference_objects(filename):
    name_to_object = {}
    for line in open(filename):
        windows_formatted_lines = line.split('\r')
        for wl in windows_formatted_lines:
            name, n20, sequence, pam = wl.split('\t')
            name_to_object[name] = Reference(name, n20, sequence, pam)
    return name_to_object

def count_indels(reference_name_to_reads, reference_name_to_template):
    for reference_name in reference_name_to_reads.keys():
        reads = reference_name_to_reads[reference_name]
        template = reference_name_to_template[reference_name]
        all_indels = get_indels_from_reads(reads, template)
        #print_and_tally_indels(all_indels, reads, template)

        plot_all_indels(all_indels, template)

def plot_all_indels(all_indels, template):

    insertions = filter(lambda x: x.is_insertion, all_indels)
    deletions = filter(lambda x: x.is_deletion, all_indels)

    insertion_dist, insertion_length, insertion_count = get_all_points_for_indels(insertions)
    deletion_dist, deletion_length, deletion_count = get_all_points_for_indels(deletions)

    # plotting
    fig = plt.figure()
    ax = fig.add_subplot(222, projection='3d')
    ax2 = fig.add_subplot(221)

    ax.scatter(insertion_dist, insertion_length, insertion_count, c='r')
    ax.scatter(deletion_dist, deletion_length, deletion_count, c='b')
    ax.set_xlabel('Distance from Cutsite')
    ax.set_ylabel('Indel Length')
    ax.set_zlabel('Number of Indels')

    insertion_distance, insertion_count = get_distance_counts_for_indels(insertions)
    deletion_distance, deletion_count = get_distance_counts_for_indels(deletions)
    ax2.scatter(insertion_distance, insertion_count, c='r')
    ax2.scatter(deletion_distance, deletion_count, c='b')
    ax2.set_xlabel('Distance from Cutsite')
    ax2.set_ylabel('Number of Indels')

    plt.savefig(template.name + ".png")

def get_all_points_for_indels(all_indels):
    distance_length_to_count = {}
    for i in all_indels:
        distance_length = (i.distance_to_cutsite, i.length)
        if distance_length in distance_length_to_count:
             distance_length_to_count[distance_length] += 1
        else:
            distance_length_to_count[distance_length] = 1

    distances = []
    lengths = []
    counts = []
    for distance_length, count in distance_length_to_count.iteritems():
        distance, length = distance_length
        distances.append(distance)
        lengths.append(length)
        counts.append(count)

    return distances, lengths, counts

def get_distance_counts_for_indels(all_indels):
    distance_to_count = {}
    for i in all_indels:
        if i.distance_to_cutsite in distance_to_count:
             distance_to_count[i.distance_to_cutsite] += 1
        else:
            distance_to_count[i.distance_to_cutsite] = 1

    distances = distance_to_count.keys()
    counts = [distance_to_count[k] for k in distances]

    return distances, counts

def print_and_tally_indels(all_indels, reads, template):
    print "name:", template.name
    print "cutsite index:", template.cutsite_index()
    print "pam:", template.pam
    print "length of reference sequence:", len(template.sequence)
    print "distance to cutsite, number of indels:"
    indel_distances_count = tally_distance_count(all_indels)
    for key in sorted(indel_distances_count.keys()):
        print key, ",", indel_distances_count[key]
    print "number of indels per read:", len(all_indels) * 1.0/len(reads)
    print "number of reads:", len(reads)
    print "----"

def tally_distance_count(indels):
    distance_to_count = {}
    for i in indels:
        if i.distance_to_cutsite in distance_to_count:
            distance_to_count[i.distance_to_cutsite] += 1
        else:
            distance_to_count[i.distance_to_cutsite] = 1
    return distance_to_count

def run():
    sam_file = "output.sam"
    template_file = "refseq_guide.txt"
    reference_name_to_reads = parse_sam(sam_file)
    reference_name_to_template = get_reference_objects(template_file)

    reads = parse_sam(sam_file)
    references = parse_references(template_file)

    count_indels(reference_name_to_reads, reference_name_to_template)

run()
