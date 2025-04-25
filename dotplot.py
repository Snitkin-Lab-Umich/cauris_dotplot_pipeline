import os
import subprocess
import argparse


def check_query(input_path):
    query_fasta_list = []
    if os.path.isdir(input_path):
        for f in os.listdir(input_path):
            if check_suffix(f):
                query_fasta_list.append(input_path + f)
    elif os.path.isfile(input_path) and check_suffix(input_path):
        query_fasta_list.append(input_path)
    else:
        print('Unrecognized path for ')
        quit(1)
    return(query_fasta_list)

def check_suffix(path,suffix_list = ['.fasta','.fa','.fna']):
    if any([path.endswith(x) for x in suffix_list]):
        return(True)
    else:
        return(False)

def run_nucmer(query_list,subject,output_dir,debug):
    output_file_list = []
    subject_name = subject.split('/')[-1].split('.fa')[0]
    with open(debug,'a') as debug_log:
        for query in query_list:
            query_name = query.split('/')[-1].split('.fa')[0]
            file_prefix = f'{query_name}_to_{subject_name}'
            # create the alignment file (A_to_B.delta)
            command1 = ['nucmer','-p',file_prefix,subject,query]
            subprocess.run(command1)
            _ = debug_log.write(' '.join(command1)+'\n')
            # create the coordinate file (A_to_B.coord)
            with open(file_prefix + '.coord','w') as fh_coord:
                command2 = ['show-coords','-r','-c','-l','-T',file_prefix + '.delta']
                subprocess.run(command2,stdout = fh_coord)
                _ = debug_log.write('_'.join(command1)+'\n')
            # move both files to the results directory
            subprocess.run(['mv',file_prefix + '.delta',output_dir])
            subprocess.run(['mv',file_prefix + '.coord',output_dir])
            # append the final output to the list
            output_file_list.append(file_prefix + '.coord')
    return(output_file_list)

def count_contig_len(input_file,output_file):
    with open(input_file,'r') as fh:
        t = None
        data_list = []
        for line in fh:
            line = line.strip()
            if line[0] == '>':
                if t is not None:
                    data_list.append(t)
                t = [line.split('>')[1],0]
            else:
                t[1]+=len(line)
    data_list.append(t)
    with open(output_file,'w') as fh:
        _ = fh.write('contig_name\tcontig_length\n')
        for el in data_list:
            #print('\t'.join([str(x) for x in el])+'\n')
            _ = fh.write('\t'.join([str(x) for x in el])+'\n')

def make_plots(nucmer_dir,nucmer_files,contig_data_dir,highlight_data,output_dir):
    for filename in nucmer_files:
        query_name,subject_name = filename.split('.coord')[0].split('_to_')
        query_contig_data = contig_data_dir + query_name + '_contig_data.csv'
        subject_contig_data = contig_data_dir + subject_name + '_contig_data.csv'
        output_file = output_dir + f'{query_name}_to_{subject_name}.pdf'
        command = ['Rscript','make_plots.R',nucmer_dir + filename,subject_contig_data,query_contig_data,highlight_data,output_file]
        subprocess.run(command)

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--query','-q',type=str,
        help='''Provide a query to perform the alignments to. This can be either a single file or a directory of files.
        All subject files should be in the fasta format (.fasta or .fa).''',
        default=None,required=True
        )
    parser.add_argument(
        '--subject','-s',type=str,
        help='''Provide a subject sequence in the nucleotide fasta format (.fasta, .fa, or .fna). This should be a single file.''',
        default=None,required=True
        )
    parser.add_argument(
        '--name','-n',type=str,help='''(Optional) Provide a name for the run. This name will be used for all outputs.''',
        default='new_search'
        )
    parser.add_argument(
        '--highlight','-hl',type=str,help='''(Optional) Provide a file containing regions to highlight on the plot. This needs to 
        match the format of the provided highlight_data.tsv file.''',
        default='NA'
        )
    # parser.add_argument(
    #     '--container','-c',type=str,help='''Provide a singularity container for MUMmer.
    #     The specified container will be used, instead of assuming you have MUMmer installed locally.''',
    #     default=None
    #     )
    args = parser.parse_args()
    if args.query is None or args.subject is None:
        print('Missing query or subject sequence')
        quit(1)
    # change working directory to location of script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # generate the files and directories needed for the BLAST search in this directory
    for d in ['logs/',f'temp/{args.name}/',f'results/{args.name}/',f'plots/{args.name}/']:
        if not os.path.isdir(d):
            subprocess.call(['mkdir','-p',d])
    # create log file in this directory
    with open('logs/debug_log.txt','w') as fh:
        _ = fh.write(f'debug log for {args.name}\n')
    # check that subject and query exist
    for f in [args.query,args.subject]:
        if f is not None:
            if not os.path.exists(f):
                print(f'Could not locate file or directory at {f}')
                quit(1)
    # fix paths
    #args.query,args.subject = [os.path.abspath(x) if x is not None else x for x in [args.query,args.subject]]
    # create list of query files
    query_fasta_list = check_query(args.query)
    #print(query_fasta_list)
    # determine contig lengths for the query and subject
    for fpath in query_fasta_list + [args.subject]:
        fname = fpath.split('/')[-1].split('.fa')[0]
        #print(fname)
        count_contig_len(fpath,f'temp/{args.name}/' + fname + '_contig_data.csv')
    # run nucmer, aligning each query to the subject
    coord_file_list = run_nucmer(query_list=query_fasta_list, subject=args.subject, output_dir=f'results/{args.name}/', debug='logs/debug_log.txt')
    # create the plots using the R script
    make_plots(
        nucmer_dir=f'results/{args.name}/',nucmer_files=coord_file_list,contig_data_dir=f'temp/{args.name}/',
        highlight_data=args.highlight,output_dir=f'plots/{args.name}/')


if __name__ == "__main__":
    main()







