from celery import Celery
from celery import states
from celery.exceptions import SoftTimeLimitExceeded
from flask import Flask, render_template, request, flash, redirect, url_for
from indel_mapper_lib.aws_upload import AwsUpload
from indel_mapper_lib.csv_upload import CsvUpload
from indel_mapper_lib.csv_upload import NullCsvUpload
from indel_mapper_lib.csv_writer import CsvWriter
from indel_mapper_lib.json_writer import JsonWriter
from indel_mapper_lib.presigned_post import PresignedPost
from indel_mapper_lib.presenter import Presenter
from indel_mapper_lib.reference_parser import ReferenceParser
from indel_mapper_lib.reference_presenter import ReferencePresenter
from indel_mapper_lib.sam_parser import SamParser
import binascii
import csv
import json
import os
import pysam
import urllib.request

LOCAL_REDIS_URL = "redis://localhost:6379/0"

# Flask

app = Flask(__name__)
app.secret_key = os.environ["SECRET_FLASK_KEY"]
app.config["UPLOAD_FOLDER"] = "/tmp"

if "S3_ACCESS_KEY" in os.environ and "S3_SECRET_KEY" in os.environ:
    app.s3_access_key = os.environ["S3_ACCESS_KEY"]
    app.s3_secret_key = os.environ["S3_SECRET_KEY"]
    app.s3_is_configured = True
else:
    app.s3_is_configured = False

if "REDIS_URL" in os.environ:
    app.config["CELERY_BROKER_URL"] = os.environ["REDIS_URL"]
    app.config["CELERY_RESULT_BACKEND"] = os.environ["REDIS_URL"]
else:
    app.config["CELERY_BROKER_URL"] = LOCAL_REDIS_URL
    app.config["CELERY_RESULT_BACKEND"] = LOCAL_REDIS_URL

# Celery

celery = Celery(app.name,
                broker=app.config["CELERY_BROKER_URL"],
                task_track_started=True,
                task_soft_time_limit=350)
celery.conf.update(app.config)

@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")

@celery.task(bind=True)
def compute_indels_near_cutsite(self, csv_path, sam_path):
    try:
        if app.s3_is_configured:
            tmp_csv_path, _ = urllib.request.urlretrieve(csv_path)
            tmp_sam_path, _ = urllib.request.urlretrieve(sam_path)
        else:
            tmp_csv_path = csv_path
            tmp_sam_path = sam_path

        reference_name_to_reads = SamParser(pysam.AlignmentFile(tmp_sam_path, "rb")).reference_name_to_reads_dict()
        references = ReferenceParser(csv.reader(open(tmp_csv_path)), reference_name_to_reads).references()
        presenter_results = Presenter([reference for reference in references if reference.is_valid])
        results = presenter_results.present()

        # Write a temporary results CSV file.
        tmp_results_csv_path = _random_tempfile_path()
        CsvWriter(results).write_to(tmp_results_csv_path)
        # Write a temporary json to store display information.
        tmp_json_path = _random_tempfile_path()
        JsonWriter(results).write_to(tmp_json_path)

        os.remove(tmp_sam_path)
        os.remove(tmp_csv_path)

        if app.s3_is_configured:
            uploaded_results = _upload_results(tmp_results_csv_path)
            uploaded_json = _upload_file(tmp_json_path, _random_filename())

            os.remove(tmp_results_csv_path)
            os.remove(tmp_json_path)

            return {"csv": uploaded_results.url, "json": uploaded_json.url}
        else:
            return {"csv": tmp_results_csv_path, "json": tmp_json_path}
    except SoftTimeLimitExceeded:
        os.remove(tmp_sam_path)
        os.remove(tmp_csv_path)

        if app.s3_is_configured:
            os.remove(tmp_results_csv_path)
            os.remove(tmp_json_path)

        return {"csv": "", "json": ""}

def _random_filename():
    return binascii.hexlify(os.urandom(16)).decode("utf-8")

def _random_tempfile_path():
    return os.path.join(app.config["UPLOAD_FOLDER"], _random_filename())

def _upload_file(local_filename, remote_filename):
    return AwsUpload(app.s3_access_key,
                    app.s3_secret_key,
                    open(local_filename, "rb"),
                    remote_filename)

def _upload_results(results_csv_filename):
    if not app.s3_is_configured:
        print("AWS isn't configured, can't upload CSV to S3.")
        return NullCsvUpload()

    # Upload it to S3.
    upload = CsvUpload(app.s3_access_key, app.s3_secret_key,
                       open(results_csv_filename, "rb"))

    return upload

@app.route('/status/<task_id>')
def taskstatus(task_id):
    task = compute_indels_near_cutsite.AsyncResult(task_id)
    if task.state == states.SUCCESS:
        task_csv = task.result["csv"]
        task_json = task.result["json"]

        if task_csv == "" or task_json == "":
            flash("Error processing.")
            return render_template("status.html", results=[], upload="", processing=False)

        if app.s3_is_configured:
            tmp_json_path, _ = urllib.request.urlretrieve(task_json)
        else:
            tmp_json_path = task_json

        json_file = open(tmp_json_path, "r")
        reference_presenter_results = json.loads(json_file.read())

        os.remove(tmp_json_path)

        return render_template("status.html", results=reference_presenter_results, upload=task_csv, processing=False)
    elif task.state == states.FAILURE:
        flash("Error processing.")
        return render_template("status.html", results=[], upload="", processing=False)
    else:
        return render_template("status.html", results=[], upload="", processing=True)

@app.route('/sign_s3/', methods=['GET'])
def sign_s3():
    if not app.s3_is_configured:
        return json.dumps({"message": "AWS not configured."}), 500

    presigned_post = PresignedPost(app.s3_access_key,
                                   app.s3_secret_key,
                                   _random_filename(),
                                   request.args.get('file_type'))

    if presigned_post.succeeded:
        return json.dumps({
            "data": presigned_post.presigned_data,
            "url": presigned_post.url
        })
    else:
        return json.dumps({"message": "Could not get signed request."}), 500

def isNoneOrEmpty(url_from_request_parameter):
    return url_from_request_parameter is None or url_from_request_parameter == ""

@app.route('/process_files/', methods=['POST'])
def process_results():
    reference_url = request.form.get('reference-hidden')
    alignment_url = request.form.get('alignment-hidden')

    if not isNoneOrEmpty(reference_url) and not isNoneOrEmpty(alignment_url):
        results = compute_indels_near_cutsite.apply_async(args=[reference_url, alignment_url])
        # Display the status page containing the results
        print(results.id)
        return redirect(url_for('taskstatus', task_id=results.id))
    else:
        flash("Missing files.")
        return redirect(url_for('index'))

if __name__ == "__main__":
    app.run()
