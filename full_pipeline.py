import argparse

parser = argparse.ArgumentParser(description='Process raw R2C2 read data into csv files containing payloads + data')

parser.add_argument('read_dir', type=str, help='Input folder containing files with raw read data. All fastq files in the given directory will be processed')
parser.add_argument('output', type=str, help='Directory path for all output data. Should point to an empty directory')

args = parser.parse_args()

read_dir = args.read_dir
output_dir = args.output

import os

def get_files(read_dir, suffix):
    if os.path.exists(read_dir):
        file_list = os.listdir(read_dir) # Get list of all files in directory
        file_list = [e in file_list if e.endswith(suffix)] # Choose files that end with the proper suffix
        return file_list
    else:
        return -1

def pre_cmd(cp3oa_pre, read_file, out_path, cut_qual, cut_len, splint_seq, config):
    cmd = 'python {} -i {} -o {} -q {} -l {} -s {} -c {}'.format(c3poa_pre, read_file, out_path, cut_qual, cut_len, splint_seq, config)
    os.system(cmd)

def preprocess(output_dir, c3poa_pre, read_dir, cut_qual, cut_len, splint_seq, config):
    # Set output directory
    p_out = output_dir + '/pre'
    if not os.path.exists(p_out):
        os.makedirs(p_out)
    # Get input files
    file_list = [e in os.listdir(read_dir) if e.endswith('.fastq')]
    # Run preprocess on each file in the directory
    for fname in file_list:
        read_file = read_dir + '/' + fname
        pre_cmd(c3poa_pre, read_file, p_out, cut_qual, cut_len, splint_seq, config) 
    return p_out

def c3poa_cmd(c3poa, read_file, consensus_out, temp, mat, l, d, config):
    cmd = 'python {} -r {} -p {} -m {} -l {} -d {} -c {} -o {}'.format(c3poa, read_file, temp, mat, l, d, config, consensus_out)
    os.system(cmd) # Run command

def cproc_sdir(sd_name, d_label, mat, l, d, config):
    label = sd_name[-2:] # Get directory label (which splint is it?)
    consensus_out = c_out_file.format(label) # Determine consensus out file
    f_name = 'R2C2_raw_reads.fastq' # Name of input file --FLAG: MAGIC NUMBER--
    if f_name in os.listdir(sd_name): # If input file in splint directory:
        read_file = sd_name + '/' + f_name # Set name of input file (with directory filepath)
        temp = c_temp.format(d_label, label) # Set temp file name
        if not os.path.exists(temp): # Make temp directory if it doesn't exist
            os.makedirs(temp)
        c3poa_cmd(process, read_file, consensus_out, temp, mat, l, d, config)

def cproc_dir(dir_name, c_out_file, c_temp, mat, l, d, config):
    d_label = dir_name.split('/')[-1] # Get directory number
    # Get splint directories for preproc directory
    splint_dirs = [(dir_name + '/' + name) for name in os.listdir(dir_name) if os.path.isdir(dir_name + '/' + name)]
    for sd_name in splint_dirs: # Iterate through splint directories
        cproc_sdir(sd_name, d_label, mat, l, d, config)

def run_c3poa(output_dir, process, pre_path, mat, l, d, config): 
    # Set c3poa output files
    c_out = output_dir + '/consensus' # C3POa output
    c_out_file = c_out + '/R2C2_consensus_g{}.fasta'
    c_temp = c_out + '/temp_{}_{}'
    # Get preprocessing directory paths
    preprocess_dirs = [(pre_path + '/' + name) for name in os.listdir(pre_path) if os.path.isdir(pre_path + '/' + name)]
    for d in preprocess_dirs: # Iterate through preprocessing directories
        cproc_dir(d, c_out_file, t_temp, mat, l, d, config)
    return c_out, c_out_file

# pipe files from /consensus into splint_remove.py | pipe splint_remove.py into /payload folder

def splint_cmd(splint_remove, sr_input, s_path, sr_temp, sr_out):
    cmd = 'python {} {} {} {} {}'.format(splint_remove, sr_input, s_path, sr_temp, sr_out)
    os.system(cmd)

def splint_remove(output_dir, splint_remove, splint_path, c_out_file):
    s_out = output_dir + '/payload' # splint_remove output
    s_outf = s_out + '/payload_g{}.fasta'
    s_outc = s_out + '/payload_g{}.csv'
    s_temp = s_out + '/temp_{}'
    for i in range(0, 25): # --FLAG: Magic number? Could generalize. Set aside for later.
        label = str(i + 1)
        if int(label) < 10:
            label = '0' + label
        s_path = splint_path.format(label)
        sr_input = c_out_file.format(label)
        sr_temp = s_temp.format(label)
        sr_outf = s_outf.format(label)
        sr_outc = s_outc.format(label)
        if os.path.exists(sr_input):
            sr_cmd(splint_remove, sr_input, s_path, sr_temp, sr_out)

def main():
    c3poa_pre = '/ssd1/home/Cas9/software/C3POa/C3POa_preprocessing.py'
    config = '/ssd1/home/Cas9/software/C3POa/config.txt'
    cut_qual = '9'
    cut_len = '1000'
    splint_seq = '/ssd1/home/Cas9/splint_sequences/splint_all.fasta'
    pre_path = preprocess(output_dir, cp3oa_pre, read_file, cut_qual, cut_len, splint_seq, config)

    process = '/ssd1/home/Cas9/software/C3POa/C3POa.py'
    mat = '/ssd1/home/Cas9/software/C3POa/NUC.4.4.mat'
    l = '500'
    d = '110'
    c_out_file = run_c3poa(output_dir, process, pre_path, mat, l, d, config)

    splint_remove = '/ssd1/home/Cas9/splint_remove.py'
    splint_path = '/ssd1/home/Cas9/splint_sequences/splint_g{}.fasta'
    post_dir = splint_remove(splint_remove, splint_path, c_out_file)

if __name__ == '__main__':
    main()
