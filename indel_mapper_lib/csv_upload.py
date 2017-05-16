import datetime
import os
import tinys3


class CsvUpload(object):
    BUCKET_NAME = "indel-mapper"
    S3_ACCESS_KEY = os.environ["S3_ACCESS_KEY"]
    S3_SECRET_KEY = os.environ["S3_SECRET_KEY"]

    def __init__(self, file):
        self._try_uploading(file)

    def _try_uploading(self, file):
        file_name = self._csv_file_name()
        try:
            conn = tinys3.Connection(self.S3_ACCESS_KEY, self.S3_SECRET_KEY, tls=True)
            conn.upload(file_name, file, self.BUCKET_NAME)
            self.succeeded = True
            self.url = self._url_for(file_name)
        except:
            self.succeeded = False
            self.url = None

    def _url_for(self, file_name):
        return "https://{}.s3.amazonaws.com/{}".format(self.BUCKET_NAME, file_name)

    def _csv_file_name(self):
        now = datetime.datetime.now()
        return now.strftime("results-%Y-%m-%d-%H-%M-%S-%f.csv")
