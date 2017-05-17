import boto3
import datetime
import os


class CsvUpload(object):
    BUCKET_NAME = "indel-mapper"
    S3_ACCESS_KEY = os.environ["S3_ACCESS_KEY"]
    S3_SECRET_KEY = os.environ["S3_SECRET_KEY"]

    def __init__(self, file):
        self._try_uploading(file)

    def _try_uploading(self, file):
        filename = self._csv_filename()
        try:
            self._upload(file, filename)
            self._make_public(filename)

            self.succeeded = True
            self.url = self._url_for(filename)
        except Exception as e:
            print(e)
            self.succeeded = False
            self.url = None

    def _s3(self):
        return boto3.resource("s3",
                              aws_access_key_id=self.S3_ACCESS_KEY,
                              aws_secret_access_key=self.S3_SECRET_KEY)

    def _upload(self, file, filename):
        self._s3().Bucket(self.BUCKET_NAME).upload_fileobj(file, filename)

    def _make_public(self, filename):
        object_acl = self._s3().ObjectAcl(self.BUCKET_NAME, filename)
        object_acl.put(ACL="public-read")

    def _url_for(self, filename):
        return "https://{}.s3.amazonaws.com/{}".format(self.BUCKET_NAME, filename)

    def _csv_filename(self):
        now = datetime.datetime.now()
        return now.strftime("results-%Y-%m-%d-%H-%M-%S-%f.csv")


class NullCsvUpload(object):
    def __init__(self):
        self.succeeded = False
        self.url = None
