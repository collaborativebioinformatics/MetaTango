import subprocess
import pathlib
import os
import pandas as pd

class SylphWrapper:

    def __init__(self):
        self.sylph = "sylph"
        self.sylph_to_tax = "sylph-tax"
        self.gtdb_data_info = {'path_val' : None, 'name' : None, 'c' : None}
        self.profile_tsv = None

    def run_cmd(self, cmd, capture_output=False):
        result = subprocess.run(cmd, check=True, capture_output=capture_output, text=True)
        return result.stdout if capture_output else None
    
    def set_gtdb_database_path(self,gtdb_path_name,gtdb_name,gtdb_c):
        # TODO: function to update inbuilt path value using string
        """set parameters to existing path"""
        self.gtdb_data_info['path_val'] = gtdb_path_name
        self.gtdb_data_info['name'] = gtdb_name
        self.gtdb_data_info['c'] = gtdb_c

    def get_gtdb_database(self,name : str = 'r214',c_val : int = 200, wget_save_dir = './data_files/'):
        """
        download gtdb database into local directory
        """
        
        if name not in ['r214', 'r226', 'r220']:
            raise ValueError()
        if c_val not in [200,1000]:
            raise ValueError()
        self.gtdb_data_info['name'] = name
        self.gtdb_data_info['c'] = c_val
        
        base_url = "http://faust.compbio.cs.cmu.edu/sylph-stuff/"
        if name in ['r226', 'r220']:
            file_name = f'gtdb-{name}-c{c_val}-dbv1.syldb'
        if name in ['r214']:
            file_name = f'v0.3-c{c_val}-gtdb-{name}.syldb'

        cmd_lst = ["wget", "-P", wget_save_dir, f"{base_url}{file_name}"]
        os.makedirs(wget_save_dir,exist_ok=True)

        self.gtdb_data_info['path_val'] = os.path.join(wget_save_dir,file_name)

        # TODO: RUN ONLY IF PATH DOESN'T EXIST
        self.run_cmd(cmd_lst)

    def sketch_single_end(self, read_fastq: str, read_fasta: str = None):
        """runs the sketch command on single-end reads"""
        # TODO: FIGURE OUT WHERE SKETCH SAVES THE .SYL FILES AND RETURN THAT PATH        
        # i think this needs the use of '-S' ?
        # TODO: FIGURE OUT WHAT TO DO IF READ_FASTA IS MULTIPLE FILES

        if read_fasta is not None:
            # sylph sketch -r reads1.fq reads2.fa
            self.run_cmd([self.sylph, "sketch", "-r", read_fastq, read_fasta])
        else:
            # sylph sketch reads1.fq 
            self.run_cmd([self.sylph, "sketch", read_fastq]) 

    def profile(self,sketched_fastq_path : str, result_tsv_path : str, num_threads : int = 4):
        """profiles the sketched fastqs by comparing them to the gtdb database downloaded / set"""

        # TODO - CAN PROBABLY READ THE SKETCHED_FASTQ_PATH FROM FN SKETCH_SINGLE END 
        # INSTEAD OF PASSING IT AS AN ARG (?)        

        self.run_cmd([self.sylph, "profile", 
                        self.gtdb_data_info['path_val'], 
                        sketched_fastq_path, 
                        '-t', num_threads,
                        '-o', result_tsv_path])
        
        self.profile_tsv = result_tsv_path
    
    def get_taxonomy(self,tax_results_name):
        """
        sylph profile only gives genomic outputs that need to be integrated
        with the taxonomy from GTDB
        outputs a new file called prefix_MYSAMPLENAME.sylphmpa for each sample in the results file
        """

        # ASSUMING THAT THIS HAS BEEN DONE BEFORE:
        # mkdir taxonomy_file_folder
        # sylph-tax download --download-to taxonomy_file_folder
        
        gtdb_name = self.gtdb_data_info['name']
        gtdb_for_tax = f'GTDB_{gtdb_name}'

        self.run_cmd([
            self.sylph_to_tax, 
            "taxprof",
            self.profile_tsv,
            "-t",
            gtdb_for_tax,
            "-o",
            tax_results_name
        ])

        self.taxonomy_tsv = tax_results_name

    def get_gtdb_strain(self,symph_file):
        """get set of the gtdb strains from the output of sylph-tax"""
        
        symph_tax_df = pd.read_csv(symph_file,sep='\t', comment='#')
        clade_lst = list(symph_tax_df['clade_name'].values)
        # Split by '|' and take the part starting with 't__'
        find_strain_from_t = lambda s : [x for x in s.split('|') if x.startswith('t__')][0].replace('t__', '')        
        strain_lst = [find_strain_from_t(v) for v in clade_lst]
        strain_set = set(strain_lst)       
        return strain_set

    def convert_taxonomy_gtdb_to_ncbi(self):
        pass

class ReferenceDownloadWrapper():
    pass

if __name__ == "__main__":

    wrapper = SylphWrapper()
    wrapper.get_gtdb_database(name="r214", c_val=200)
