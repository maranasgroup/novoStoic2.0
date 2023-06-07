import sys
import os
sys.path.append('../')

current_dir = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(
    current_dir, '../data/', 'optstoic_db_v3'))
