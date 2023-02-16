"""
python3 Monitor_for_dozor_real_time_v1.py /gpfs/cfel/group/cxi/scratch/2020/EXFEL-2019-Schmidt-Mar-p002450/scratch/galchenm/scripts_for_work/smart_chip_screening/dozor/Lyso4/grid_00
"""

import os

import sys
import h5py as h5
import numpy as np
import glob
import subprocess
import time
import re
import datetime


MIN_NUM_PEAKS = 30

hdf5path = '/entry/data/data'
MAINPATH_TO_CC = '/gpfs/cfel/group/cxi/scratch/2020/EXFEL-2019-Schmidt-Mar-p002450/scratch/galchenm/scripts_for_work/Crystal_Control/CrystalControlMaxwell'

def parsing_dozor(dozorr_file):
    
    
    NumOFline = re.search(r'\d+\.out',os.path.basename(dozorr_file)).group().split('.')[0]
    number_for_hdf5_file = NumOFline.zfill(6)
    output = os.path.dirname(dozorr_file)
    
    output_filename = os.path.join(output, NumOFline+'.lst')
    
    
    with open(dozorr_file, 'r') as f:
        start = 0
        for line in f:
             
            if '+ image_data_filename   =' in line:
                
                image_filename = line.split(' = ')[1].replace('<','').replace('>','').strip()
                image_data_filename = image_filename.replace('master',f'data_{number_for_hdf5_file}')
                h5file = h5.File(image_data_filename, 'r')
                num = h5file[hdf5path].shape[0]
                h5file.close()
                
                continue
            elif 'image | num.of   Score    Resolution' in line:
                start = 1
                continue
            elif start == 1:
                
                start = 2
                continue
            elif start != 2:
                continue
                
            elif '------------------------------------' in line and start == 2:
                break
            
            tmp = re.search(r'\d+ \|[^\S\n\t]+\d+', line).group()
            number_of_pattern, num_peaks = tmp.split('|')
            number_of_pattern = int(number_of_pattern.strip())
            num_peaks = int(num_peaks.strip())
            
            if num_peaks > MIN_NUM_PEAKS:
                if os.path.exists(output_filename):
                    append_write = 'a' # append if already exists
                else:
                    append_write = 'w' # make a new file if not
                
                imageNumber = number_of_pattern - ((int(NumOFline)-1)*num)
                out = open(output_filename, append_write)
                
                out.write(f'{imageNumber},{int(NumOFline)-1}\n')
                out.close()



def is_folder_updated(path, threshold_in_minutes = 0.5):
    last_forder_update_timestamp = os.stat(path).st_mtime

    last_folder_update_datetime = datetime.datetime.fromtimestamp(last_forder_update_timestamp)
    current_datetime = datetime.datetime.now()

    difference = current_datetime - last_folder_update_datetime
    
    return (difference.total_seconds()/60) < threshold_in_minutes


def Monitor(current_folder_with_dozor_results):
    
    list_of_dozorr_files = []
    after = len(list_of_dozorr_files)
    
    list_of_processed_dozorr_files = []
    
    #wait when the files start appearing
    while(after == 0):
        list_of_dozorr_files =  glob.glob(f'{current_folder_with_dozor_results}/dozor/dozorr*out')
        after = len(list_of_dozorr_files)
        
    while(1):
        list(map(parsing_dozor, list_of_dozorr_files))
        list_of_processed_dozorr_files += [dozorr_file  for dozorr_file in list_of_dozorr_files if dozorr_file not in list_of_processed_dozorr_files]
        
        time.sleep(2.5)
        if not is_folder_updated(current_folder_with_dozor_results):
            list_of_dozorr_files =  glob.glob(f'{current_folder_with_dozor_results}/dozor/dozorr*out')
            
            if sorted(list_of_processed_dozorr_files) == sorted(list_of_dozorr_files):
                break
    
if __name__ == '__main__':
    current_folder_with_dozor_results = sys.argv[1]
    
    while 1:
        if os.path.exists(current_folder_with_dozor_results):
            print('appeared')
            break
    
    Monitor(current_folder_with_dozor_results)
    
    files_to_cat = glob.glob(f'{current_folder_with_dozor_results}/dozor/*lst')
    
    command_line = "cat " + " ".join(files_to_cat) + f" > {MAINPATH_TO_CC}/with_duplicates-joined.lst"
        
    os.system(command_line)
    
    command_line = f'sort -u {MAINPATH_TO_CC}/with_duplicates-joined.lst > {MAINPATH_TO_CC}/joined.lst'
    os.system(command_line)