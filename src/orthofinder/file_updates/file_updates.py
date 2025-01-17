from . import hogs
import os

def update_output_files(input_dir, output_dir):
    hog_files = os.listdir(input_dir)
    n0_file = [
        os.path.join(input_dir, file) 
        for file in hog_files
        if file == "N0.tsv"
    ][0]
    name_dictionary, new_og_dict, header_str = hogs.update_hogs(n0_file)
    hogs.convert_hog2og_n0(new_og_dict, header_str, os.path.join(input_dir, "N0.tsv"))
    hogs.update_all_hog_files(input_dir, output_dir, name_dictionary)
    
