import os
import pyinotify
import sys
import h5py as h5
import numpy as np
import glob

import time

def calling_peakfinder8(file):
    global output_dir
    
    job_file = os.path.join(output_dir, os.path.basename(file).split('.')[0] + '.sh')
    print('process ', os.path.basename(file), '\n')
    with open(job_file, 'w+') as fh:
        fh.writelines("#!/bin/sh\n")
        fh.writelines("#SBATCH --job=%s.sh\n" % file.split('.')[0])
        fh.writelines("#SBATCH --partition=upex\n")
        fh.writelines("#SBATCH --time=12:00:00\n")
        fh.writelines("#SBATCH --nodes=1\n")
        fh.writelines("#SBATCH --mem=500000\n")
        fh.writelines("module load anaconda3/5.2\n")
        command = 'python3 pf8_coordinates_generating.py {} {}'.format(file, output_dir)
        fh.writelines(command)
    time.sleep(2.5)
    os.system("sbatch %s" % job_file)


def Monitor(path):
    
    list_of_files = glob.glob(f'{path}/*h5')
    processed = []
    after = len(list_of_files)
    while(after == 0):
        list_of_files =  glob.glob(f'{path}/*h5')
        after = len(list_of_files)
        
    
    for file in list_of_files:
        if file not in processed and 'data' in os.path.basename(file):
            processed.append(file)
            
            calling_peakfinder8(file)
    time.sleep(1)
    list_of_files = glob.glob(f'{path}/*h5')

    
if __name__ == '__main__':
    path = sys.argv[1]
    output_dir = sys.argv[2]
    while 1:
        
        if os.path.exists(path):
            print('appeared')
            
            break
    Monitor(path)
