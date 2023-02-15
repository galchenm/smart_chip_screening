import os
import pyinotify
import sys
import h5py as h5
import numpy as np
import glob
import datetime
import time

def calling_peakfinder8(file):
    global output_folder_after_pf8
    
    job_file = os.path.join(output_folder_after_pf8, os.path.basename(file).split('.')[0] + '.sh')
    print('process ', os.path.basename(file), '\n')
    with open(job_file, 'w+') as fh:
        fh.writelines("#!/bin/sh\n")
        fh.writelines("#SBATCH --job=%s.sh\n" % file.split('.')[0])
        fh.writelines("#SBATCH --partition=upex\n")
        fh.writelines("#SBATCH --time=12:00:00\n")
        fh.writelines("#SBATCH --nodes=1\n")
        fh.writelines("#SBATCH --mem=500000\n")
        fh.writelines("module load anaconda3/5.2\n")
        command = 'python3 pf8_coordinates_generating.py {} {}'.format(file, output_folder_after_pf8)
        fh.writelines(command)
    time.sleep(2.5)
    os.system("sbatch %s" % job_file)

def is_folder_updated(path, threshold_in_minutes = 2):
    last_forder_update_timestamp = os.stat(path).st_mtime

    last_folder_update_datetime = datetime.datetime.fromtimestamp(last_forder_update_timestamp)
    current_datetime = datetime.datetime.now()

    difference = current_datetime - last_folder_update_datetime
    
    return (difference.total_seconds()/60) < threshold_in_minutes

def Monitor(raw_folder):
    
    list_of_raw_data = []
    processed_by_pf8 = []
    after = len(list_of_raw_data)
    
    #wait when the files start appearing
    while(after == 0):
        list_of_raw_data =  glob.glob(f'{raw_folder}/*h5')
        after = len(list_of_raw_data)
        
    while(1):
        for file in list_of_raw_data:
            if file not in processed_by_pf8 and 'data' in os.path.basename(file):
                processed_by_pf8.append(file)
                calling_peakfinder8(file)
                
        time.sleep(2.5)
        list_of_raw_data = glob.glob(f'{raw_folder}/*h5')
        if sorted(processed_by_pf8) == sorted(list_of_raw_data) and not is_folder_updated(raw_folder):
            break

    
if __name__ == '__main__':
    raw_folder = sys.argv[1]
    output_folder_after_pf8 = sys.argv[2]
    
    if not os.path.exists(output_folder_after_pf8):
        os.mkdir(output_folder_after_pf8)
    
    while 1:    
        if os.path.exists(raw_folder):
            print('appeared')
            break
    Monitor(raw_folder)
