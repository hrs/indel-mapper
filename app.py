from flask import Flask, render_template, request, flash
from flask_uploads import (UploadSet, configure_uploads, UploadNotAllowed)
from indel_mapper_lib import SamParser
from indel_mapper_lib import ReferenceParser
from indel_mapper_lib import Presenter
import pysam

# Flask

app = Flask(__name__)
app.config["UPLOADS_DEFAULT_DEST"] = "storage"
app.secret_key = ('\xa3\xb6\x15\xe3E\xc4\x8c\xbaT\x14\xd1:'
              '\xafc\x9c|.\xc0H\x8d\xf2\xe5\xbd\xd5')

# uploads

uploaded_alignment_files = UploadSet('alignments', ('sam'))
uploaded_reference_files = UploadSet('references', ('txt'))

configure_uploads(app, (uploaded_alignment_files, uploaded_reference_files))

@app.route("/", methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        reference = request.files.get('reference')
        alignment = request.files.get('alignment')
        if not (alignment and reference):
            flash("You must submit an alignment SAM file and a reference CSV file.")
        else:
            try:
                alignment_file = uploaded_alignment_files.save(alignment)
            except UploadNotAllowed:
                flash("The alignment file was not a SAM file.")
            else:
                try:
                    reference_file = uploaded_reference_files.save(reference)
                except UploadNotAllowed:
                    flash("The reference file was not a CSV file.")
                else:
                    try:
                        alignment_full_path = uploaded_alignment_files.path(alignment_file)
                        reference_full_path = uploaded_reference_files.path(reference_file)
                        results = process_data(alignment_full_path, reference_full_path)
                        return render_template("index.html", results=results)
                    except:
                        flash("Error processing.")
                        return render_template("index.html", results=['1'])
    return render_template("index.html", results=['2','3'])

def process_data(sam_file, csv_file):
    reads = SamParser(pysam.AlignmentFile(sam_file, "rb")).reads()
    references = ReferenceParser(open(csv_file), reads).references()

    filtered_references = [r for r in references if r.is_valid()]#r.name == "AC4433"]#r.is_valid()]

    reference_presenter = Presenter(filtered_references)
    return reference_presenter.present()



if __name__ == "__main__":
    app.run()
