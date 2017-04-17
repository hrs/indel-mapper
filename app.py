from flask import Flask, render_template, request, flash
from indel_mapper_lib.sam_parser import SamParser
from indel_mapper_lib.reference_parser import ReferenceParser
from indel_mapper_lib.presenter import Presenter
from io import StringIO
import csv
import os
import pysam

MAX_CONTENT_BYTES = 20 * 1024 * 1024 # 20MB

# Flask

app = Flask(__name__)
app.secret_key = os.environ["SECRET_FLASK_KEY"]
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_BYTES

@app.route("/", methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        reference_file = request.files.get('reference')
        alignment_file = request.files.get('alignment')

        if not (alignment_file and reference_file):
            flash("You must submit an alignment SAM file and a reference CSV file.")
        elif not _is_extension(reference_file, ".csv"):
            flash("The reference file must be a CSV file.")
        elif not _is_extension(alignment_file, ".sam"):
            flash("The alignment file must be a SAM file.")
        else:
            try:
                results = _compute_indels_near_cutsite(alignment_file, reference_file)
                return render_template("index.html", results=results)
            except Exception as e:
                flash("Error processing.")
                return render_template("index.html", results=[])
    return render_template("index.html", results=[])

# HTTP Error 413 Request Entity Too Large
@app.errorhandler(413)
def request_entity_too_large(error):
    flash("Uploaded file is too large.")
    return render_template("index.html", results=[]), 413

# This is a first pass naive file detection.
def _is_extension(filestorage, extension):
    filename = filestorage.filename
    return filename.endswith(extension)

def _compute_indels_near_cutsite(sam_file, csv_file):
    reference_name_to_reads = SamParser(pysam.AlignmentFile(sam_file, "rb")).reference_name_to_reads_dict()
    decoded_csv_file = csv_file.read().decode("utf-8")
    references = ReferenceParser(csv.reader(StringIO(decoded_csv_file)), reference_name_to_reads).references()

    presenter_results = Presenter([reference for reference in references if reference.is_valid])

    return presenter_results.present()


if __name__ == "__main__":
    app.run()
