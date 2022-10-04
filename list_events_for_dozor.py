import os
import pyinotify
import sys
import h5py as h5
import numpy as np
import glob
import subprocess
import time
import re



MIN_NUM_PEAKS = 30

hdf5path = '/entry/data/data'


def parsing_dozor(file):
    global dir_name
    global path
    
    NumOFline = re.search(r'\d+\.out',os.path.basename(file)).group().split('.')[0]
    number_for_hdf5_file = NumOFline.zfill(6)
    output_filename = os.path.join(dir_name, os.path.basename(os.path.dirname(path)) + '_' + os.path.basename(path)) +'.lst'
    
    print(output_filename)
      
    with open(file, 'r') as f:
        start = 0
        for line in f:
             
            if '+ image_data_filename   =' in line:
                
                image_filename = line.split(' = ')[1].replace('<','').replace('>','').strip()
                image_data_filename = image_filename.replace('master',f'data_{number_for_hdf5_file}').rstrip("\x00")
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
                
                out.write(f'{image_data_filename} //{imageNumber}\n')
                out.close()
    
def Monitor(path):
    list_of_files = glob.glob(f'{path}/dozor/dozorr*out')
    
    for file in list_of_files:
        parsing_dozor(file)



    
if __name__ == '__main__':
    path = sys.argv[1]
    dir_name = sys.argv[2]
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)    

    if os.path.exists(path):
        print('start processing')
          
        Monitor(path)
    else:
        print('Something is wrong with input path')
