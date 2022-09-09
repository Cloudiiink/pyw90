import subprocess
import os

path = os.getcwd()
log  = f'{path}/auto_w90_output.txt'
subprocess.Popen(['python', 'auto_w90_fit.py'],
                 stdin  = subprocess.DEVNULL,
                 stdout = open(log, 'w'),
                 stderr = open(log, 'w'),
                 start_new_session=True)