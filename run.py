from sam_parser import SamParser
from reference_parser import ReferenceParser
from presenter import Presenter
import pysam

def run():
    sam_file = "output.sam"
    reference_file = "refseq_guide.txt"
    reads = SamParser(pysam.AlignmentFile(sam_file, "rb")).reads()
    references = ReferenceParser(open(reference_file), reads).references()

    filtered_references = [r for r in references if r.is_valid()]#r.name == "AC4433"]#r.is_valid()]

    reference_presenter = Presenter(filtered_references)
    print(reference_presenter.present())
run()
