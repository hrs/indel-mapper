from flask import Flask, render_template, request, flash, redirect, url_for
from indel_mapper_lib.sam_parser import SamParser
from indel_mapper_lib.reference_parser import ReferenceParser
from indel_mapper_lib.presenter import Presenter
from indel_mapper_lib.csv_writer import CsvWriter
from indel_mapper_lib.csv_upload import CsvUpload
from indel_mapper_lib.csv_upload import NullCsvUpload
from indel_mapper_lib.reference_presenter import ReferencePresenter
from indel_mapper_lib.aws_upload import AwsUpload
from werkzeug.utils import secure_filename
from io import StringIO
import binascii
import csv
import os
import pysam
from celery import Celery
from celery import states
import urllib.request

MAX_CONTENT_BYTES = 45 * 1024 * 1024 # 85MB
UPLOAD_FOLDER = "/tmp"

# Flask

app = Flask(__name__)
app.secret_key = os.environ["SECRET_FLASK_KEY"]
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_BYTES
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

if "S3_ACCESS_KEY" in os.environ and "S3_SECRET_KEY" in os.environ:
    app.s3_access_key = os.environ["S3_ACCESS_KEY"]
    app.s3_secret_key = os.environ["S3_SECRET_KEY"]
    app.s3_is_configured = True
else:
    app.s3_is_configured = False

if "REDIS_URL" in os.environ:
    app.config['CELERY_BROKER_URL'] = os.environ['REDIS_URL']
    app.config['CELERY_RESULT_BACKEND'] = os.environ['REDIS_URL']
else:
    app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
    app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'

# Celery

celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)


@app.route("/", methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        reference_file = request.files.get('reference')
        alignment_file = request.files.get('alignment')

        if not (alignment_file and reference_file):
            flash("You must submit an alignment file and a reference file.")
        elif not _is_extension(reference_file, ".csv"):
            flash("The reference file must be a .csv file.")
        elif not _is_extension(alignment_file, ".sam"):
            flash("The alignment file must be a .sam file.")
        else:
            try:

                # random, unlikely to collide names

                tmp_reference_filename, tmp_alignment_filename = _get_tmp_reference_and_alignment_filenames()
                tmp_reference_path = os.path.join(app.config['UPLOAD_FOLDER'], tmp_reference_filename)
                tmp_alignment_path = os.path.join(app.config['UPLOAD_FOLDER'], tmp_alignment_filename)

                reference_file.save(tmp_reference_path)
                alignment_file.save(tmp_alignment_path)

                if app.s3_is_configured:
                    reference_upload = AwsUpload(app.s3_access_key,
                                                 app.s3_secret_key,
                                                 open(tmp_reference_path, "rb"),
                                                 tmp_reference_filename)
                    alignment_upload = AwsUpload(app.s3_access_key,
                                                 app.s3_secret_key,
                                                 open(tmp_alignment_path, "rb"),
                                                 tmp_alignment_filename)
                    results = compute_indels_near_cutsite.apply_async(args=[reference_upload.url, alignment_upload.url])
                    os.remove(tmp_reference_path)
                    os.remove(tmp_alignment_path)
                else:
                    results = compute_indels_near_cutsite.apply_async(args=[tmp_reference_path, tmp_alignment_path])
                return redirect(url_for('taskstatus', task_id=results.id))

            except Exception as e:
                print(e)
                flash("Error processing.")
                return render_template("index.html", results=[], upload=NullCsvUpload())
    return render_template("index.html", results=[], upload=NullCsvUpload())

# HTTP Error 413 Request Entity Too Large
@app.errorhandler(413)
def request_entity_too_large(error):
    flash("Uploaded file is too large.")
    return render_template("index.html", results=[], upload=NullCsvUpload()), 413

# This is a first pass naive file detection.
def _is_extension(filestorage, extension):
    filename = filestorage.filename
    return filename.endswith(extension)

def _get_tmp_reference_and_alignment_filenames():
    random_number = binascii.hexlify(os.urandom(16)).decode("utf-8")
    return "reference_{}.csv".format(random_number), "alignment_{}.sam".format(random_number)

@celery.task(bind=True)
def compute_indels_near_cutsite(self, csv_path, sam_path):
    if app.s3_is_configured:
        tmp_reference_filename, tmp_alignment_filename = _get_tmp_reference_and_alignment_filenames()
        tmp_csv_path = os.path.join(app.config['UPLOAD_FOLDER'], tmp_reference_filename)
        tmp_sam_path = os.path.join(app.config['UPLOAD_FOLDER'], tmp_alignment_filename)
        urllib.request.urlretrieve(csv_path, tmp_csv_path)
        urllib.request.urlretrieve(sam_path, tmp_sam_path)
    else:
        tmp_csv_path = csv_path
        tmp_sam_path = sam_path

    reference_name_to_reads = SamParser(pysam.AlignmentFile(tmp_sam_path, "rb")).reference_name_to_reads_dict()
    references = ReferenceParser(csv.reader(open(tmp_csv_path)), reference_name_to_reads).references()
    presenter_results = Presenter([reference for reference in references if reference.is_valid])
    results = presenter_results.present()

    json_results = [result.to_dict() for result in results]

    os.remove(tmp_sam_path)
    os.remove(tmp_csv_path)

    return json_results

def _upload(results):
    if not app.s3_is_configured:
        print("Aws isn't configured, can't upload CSV to S3.")
        return NullCsvUpload()

    # Write a temporary CSV file with an unlikely-to-collide name.
    random_string = binascii.hexlify(os.urandom(16)).decode("utf-8")
    csv_temp_filename = "/tmp/{}.csv".format(random_string)
    CsvWriter(results).write_to(csv_temp_filename)

    # Upload it to S3.
    upload = CsvUpload(app.s3_access_key, app.s3_secret_key,
                       open(csv_temp_filename, "rb"))
    os.remove(csv_temp_filename)
    return upload

@app.route('/status/<task_id>')
def taskstatus(task_id):
    task = compute_indels_near_cutsite.AsyncResult(task_id)
    if task.state == states.SUCCESS:
        reference_presenter_results = []
        for reference_presenter_dict in task.result:
            reference_presenter_results.append(ReferencePresenter(from_dict=reference_presenter_dict))
        upload = _upload(reference_presenter_results)
        return render_template("status.html", results=reference_presenter_results, upload=upload, processing=False)
    elif task.state == states.FAILURE:
        flash("Error processing.")
        return render_template("status.html", results=[], upload=NullCsvUpload(), processing=False)
    else:
        return render_template("status.html", results=[], upload=NullCsvUpload(), processing=True)

if __name__ == "__main__":
    app.run()
