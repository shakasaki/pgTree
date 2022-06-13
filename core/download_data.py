import os
from core import DATA_DIR
from pooch import __version__, check_version, HTTPDownloader, retrieve
import zipfile

url = "https://polybox.ethz.ch/index.php/s/7wXwtlmkw8V1ac2/download"
file_path = DATA_DIR + "data.zip"
url = url.format(check_version(__version__, fallback="main"))
downloader = HTTPDownloader()
downloader(url=url, output_file=file_path, pooch=None)

while not os.path.exists(DATA_DIR + os.sep + "data.zip"):
    print('File being downloaded')

file_name = os.path.abspath(file_path)  # get full path of files
zip_ref = zipfile.ZipFile(file_name)  # create zipfile object
zip_ref.extractall(DATA_DIR)  # extract file to dir
os.remove(file_path) # remove zip file