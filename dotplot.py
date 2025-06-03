import os
import subprocess
import argparse


def check_query(input_path):
    query_fasta_list = []
    if os.path.isdir(input_path):
        if not input_path.endswith('/'):
            input_path+='/'
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

def suffix_trim(fname):
    if fname.endswith('.fasta') or fname.endswith('.fa'):
        return(fname.split('.fa')[0])
    elif fname.endswith('.fna'):
        return(fname.split('.fna')[0])
    else:
        print('Incorrect file name when generating contig data')
        quit(1)


def run_nucmer(query_list,subject,output_dir,debug,rname):
    output_file_list = []
    subject_name = subject.split('/')[-1].split('.fa')[0]
    with open(debug,'a') as debug_log:
        for query in query_list:
            query_name = query.split('/')[-1].split('.fa')[0]
            #file_prefix = f'{query_name}_to_{subject_name}'
            file_prefix = rname
            # create the alignment file (A_to_B.delta)
            command1 = ['nucmer','-p',file_prefix,subject,query]
            subprocess.run(command1)
            _ = debug_log.write(' '.join(command1)+'\n')
            # create the coordinate file (A_to_B.coord)
            with open(file_prefix + '.coord','w') as fh_coord:
                command2 = ['show-coords','-r','-c','-l','-T',file_prefix + '.delta']
                subprocess.run(command2,stdout = fh_coord)
                _ = debug_log.write(' '.join(command1)+'\n')
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
                contig_name = line.split('>')[1]
                # remove anything after a space to keep contig names consistent with mummer and BLAST
                contig_name = contig_name.split(' ')[0]
                t = [contig_name,0]
            else:
                t[1]+=len(line)
    data_list.append(t)
    with open(output_file,'w') as fh:
        _ = fh.write('contig_name\tcontig_length\n')
        for el in data_list:
            _ = fh.write('\t'.join([str(x) for x in el])+'\n')

def get_nucmer_coord(nucmer_dir):
    # take a nucmer results directory and generate a list of .coord files (the same output as run_nucmer above)
    # returns and empty list if the path is not a directory
    output_file_list = []
    if os.path.isdir(nucmer_dir):
        for fname in os.listdir(nucmer_dir):
            if fname.endswith('.coord'):
                output_file_list.append(fname)
    return(output_file_list)


def make_plots(nucmer_dir,nucmer_files,contig_data_dir,highlight_data,output_dir):
    for filename in nucmer_files:
        #query_name,subject_name = filename.split('.coord')[0].split('_to_')
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
        help='''Provide a query to perform the alignments to. This can be only be a single file!
        All files should be in the fasta format (.fasta or .fa).''',
        default=None
        )
    parser.add_argument(
        '--subject','-s',type=str,
        help='''Provide a subject sequence in the nucleotide fasta format (.fasta, .fa, or .fna). This should be a single file.''',
        default=None
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
    parser.add_argument(
        '--alignments','-a',type=str,help='''(Optional) Provide a path to a previously-generated results directory. This will remake the plots using the 
        Nucmer alignments and contig data present in the directory.''',
        default=None
        )
    args = parser.parse_args()
    # these statements don't cover all possbile inputs
    if (args.query is None or args.subject is None) and (args.alignments is None):
        print('No query+subject or results directory provided')
        quit(1)
    if (args.query is not None and args.subject is not None) and (args.alignments is not None):
        print('Please provide either a query+subject or a results directory with Nucmer alignments, not both.')
        quit(1)
    # change working directory to location of script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # define all directories
    #nucmer_output_dir = f'results/{args.name}/nucmer/'
    #contig_data_dir = f'results/{args.name}/contig_data/'
    #plot_output_dir = f'plots/{args.name}/'
    #debug_log_dir = 'logs/'
    #debug_log_file = f'logs/{args.name}_debug_log.txt'
    nucmer_output_dir = f'{args.name}/nucmer/'
    contig_data_dir = f'{args.name}/contig_data/'
    plot_output_dir = f'{args.name}/plots/'
    debug_log_dir = f'{args.name}/'
    debug_log_file = f'{args.name}/{args.name}_debug_log.txt'    
    # generate the directories
    dirlist = [debug_log_dir,contig_data_dir,nucmer_output_dir,plot_output_dir]
    if args.alignments is not None:
        dirlist = [debug_log_dir,plot_output_dir]
    for d in dirlist:
        if not os.path.isdir(d):
            subprocess.call(['mkdir','-p',d])
    # create log file
    with open(debug_log_file,'w') as fh:
        _ = fh.write(f'debug log for {args.name}\n')
    # if the alignments aren't already provided, generate them
    if args.alignments is None:
        # check that subject and query exist
        for f in [args.query,args.subject]:
            if not os.path.exists(f):
                print(f'Could not locate file or directory at {f}')
                quit(1)
        # create list of query files
        query_fasta_list = check_query(args.query)
        # determine contig lengths for the query and subject
        for fpath in query_fasta_list + [args.subject]:
            #fname = fpath.split('/')[-1].split('.fa')[0]
            fname = suffix_trim(fpath)
            count_contig_len(fpath,contig_data_dir + fname + '_contig_data.csv')
        # run nucmer, aligning each query to the subject
        coord_file_list = run_nucmer(query_list=query_fasta_list, subject=args.subject, output_dir=nucmer_output_dir, debug=debug_log_file, rname=args.name)
        # create the plots using the R script
    # if the alignments are provided, make sure everything looks as expected
    else:
        if os.path.isdir(args.alignments):
            if not args.alignments.endswith('/'):
                args.alignments+='/'
            # replace paths with new versions
            nucmer_output_dir = args.alignments + 'nucmer/'
            contig_data_dir = args.alignments + 'contig_data/'
            # get the list of .coord files
            coord_file_list = get_nucmer_coord(nucmer_output_dir)
            if coord_file_list == []:
                print(f'No .coord files located in {args.alignments}')
                quit(1)
        else:
            print(f'Could not locate results directory at {args.alignments}')
            quit(1)
    # make plots using all current paths
    make_plots(
        nucmer_dir=nucmer_output_dir,nucmer_files=coord_file_list,contig_data_dir=contig_data_dir,
        highlight_data=args.highlight,output_dir=plot_output_dir,rname=args.name)
    print(f'Finished making plots!')
    print(f'Location of plots: {plot_output_dir}')
    print(f'Location of Nucmer alignments: {nucmer_output_dir}')


if __name__ == "__main__":
    main()







