__version__ = '0.0.1'

import os
PACKAGE_DIR = os.path.dirname(__file__) + os.sep
DATA_DIR = os.path.abspath(PACKAGE_DIR + '../data') + os.sep
OUTPUT_DIR = os.path.abspath(PACKAGE_DIR + '../output') + os.sep
TEST_DIR = os.path.abspath(PACKAGE_DIR + '../test') + os.sep
MODEL_DIR = os.path.abspath(PACKAGE_DIR + '../model') + os.sep

# Create folders if not existing
if not os.path.isdir(DATA_DIR):
    os.mkdir(DATA_DIR)
if not os.path.isdir(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)

