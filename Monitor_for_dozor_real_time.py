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


def parsing_dozor(dozorr_file):
    global output_folder_for_post_processing_dozorr
    #print(dozorr_file)
    NumOFline = re.search(r'\d+\.out',os.path.basename(dozorr_file)).group().split('.')[0]
    number_for_hdf5_file = NumOFline.zfill(6)
    output = os.path.join(output_folder_for_post_processing_dozorr, os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(dozorr_file)))) + '_' + os.path.basename(os.path.dirname(os.path.dirname(dozorr_file))))
    
    if not(os.path.exists(output)):
        os.mkdir(output)
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



def is_folder_updated(path, threshold_in_minutes = 2):
    last_forder_update_timestamp = os.stat(path).st_mtime

    last_folder_update_datetime = datetime.datetime.fromtimestamp(last_forder_update_timestamp)
    current_datetime = datetime.datetime.now()

    difference = current_datetime - last_folder_update_datetime
    
    return (difference.total_seconds()/60) < threshold_in_minutes


def Monitor(current_folder_with_dozor_results):
    global output_folder_for_post_processing_dozorr
    
    list_of_dozorr_files = []
    after = len(list_of_dozorr_files)
    
    list_of_processed_dozorr_files = []
    
    #wait when the files start appearing
    while(after == 0):
        list_of_dozorr_files =  glob.glob(f'{current_folder_with_dozor_results}/dozor/dozorr*out')
        after = len(list_of_dozorr_files)
        
    while(1):
        for dozorr_file in list_of_dozorr_files:
            if dozorr_file not in list_of_processed_dozorr_files:
                list_of_processed_dozorr_files.append(dozorr_file)
                parsing_dozor(dozorr_file)
                
        time.sleep(2.5)
        list_of_dozorr_files =  glob.glob(f'{current_folder_with_dozor_results}/dozor/dozorr*out')
        if sorted(list_of_processed_dozorr_files) == sorted(list_of_dozorr_files) and not is_folder_updated(current_folder_with_dozor_results):
            break
    
    output = os.path.join(output_folder_for_post_processing_dozorr, os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(list_of_dozorr_files[0])))) + '_' + os.path.basename(os.path.dirname(os.path.dirname(list_of_dozorr_files[0]))))
    return output
    
if __name__ == '__main__':
    current_folder_with_dozor_results = sys.argv[1]
    output_folder_for_post_processing_dozorr = sys.argv[2]
    
    if not os.path.exists(output_folder_for_post_processing_dozorr):
        os.mkdir(output_folder_for_post_processing_dozorr)
    
    while 1:
        if os.path.exists(current_folder_with_dozor_results):
            print('appeared')
            break
    
    folder_with_dozorr_results = Monitor(current_folder_with_dozor_results)
    files_to_cat = glob.glob(f'{folder_with_dozorr_results}/*lst')
    
    command_line = "cat " + " ".join(files_to_cat) + " > joined.lst"
        
    os.system(command_line)