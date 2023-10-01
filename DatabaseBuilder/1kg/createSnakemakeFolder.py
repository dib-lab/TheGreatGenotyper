import sys
import os
import shutil
import subprocess


def run_command(command: str):
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True, shell=True)
    except subprocess.CalledProcessError as e:
        print("Command stdout:")
        print(e.stdout)
        print("Command stderr:")
        print(e.stderr)
        raise RuntimeError(f"Command '{command}' returned non-zero exit status {e.returncode}") from e


def read_sample_ids(file_path):
    with open(file_path, 'r') as f:
        return [line.strip() for line in f.readlines()]

def read_remote_file(remote_file_path):
    with open(remote_file_path, 'r') as f:
        paths = f.readlines()
    return paths

def create_directories(base_path):
    os.makedirs(os.path.join(base_path, "results/inputs/"), exist_ok=True)

def copy_config_files(base_path):
    shutil.copy("project_config.yaml", os.path.join(base_path, "project_config.yaml"))
    shutil.copy("env.yaml", os.path.join(base_path, "env.yaml"))
    shutil.copy("coverage.tsv", os.path.join(base_path, "coverage.tsv"))
    shutil.copy("Snakefile.build_3.smk", os.path.join(base_path, "Snakefile"))
    with open(os.path.join(base_path, "config.yaml"), 'w') as f:
        outPath= os.path.join( "results/")
        tmpPath= os.path.join( "scratch/")
        f.write(f"outputFolder: {outPath}\n")
        f.write(f"tempFolder: {tmpPath}\n")
        f.write(f"kSize: 31\n")
        f.write(f"batchSize: 16\n")

def write_sample_table(base_path, sample_dict):
    with open(os.path.join(base_path, "sample_table.csv"), 'w') as f:
        f.write("sample_name,cluster\n")
        samples = [ s+ ",cluster" for s in sample_dict]
        f.write("\n".join(samples) + "\n")

def write_subsample_table(base_path, sample_dict):
    with open(os.path.join(base_path, "subsample_table.csv"), 'w') as f:
        f.write("sample_name,file\n")
        for sample, files in sample_dict.items():
            for file in files:
                f.write(f"{sample},{file}\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: python createSnakemakeFolder.py <sample_ids_file> <output_folder>")
        sys.exit(1)

    sample_ids_file = sys.argv[1]
    output_folder = sys.argv[2]

    sample_ids = read_sample_ids(sample_ids_file)
    create_directories(output_folder)
    copy_config_files(output_folder)

    sample_dict = {sample_id: [] for sample_id in sample_ids}
    
    remote_file_path = os.path.join( "remote_UCDavis_GoogleDr.files")
    paths = read_remote_file(remote_file_path)

    for path in paths:
        filename = os.path.basename(path).strip()
        sample_name = filename.split(".")[0]
        extension = filename[len(sample_name):]
        
        if extension in ['.kmer_counts.gz', '.fasta.gz', '.histo']:

            if sample_name in sample_dict:
                sample_dict[sample_name].append(path.strip())
    
    write_sample_table(output_folder , sample_dict)
    write_subsample_table(output_folder , sample_dict)
    subsample_table={}
    # Generate bash commands for touch
    for sample, files in sample_dict.items():
        subsample_table[sample] = []
        for file_path in files:
            filename = os.path.basename(file_path).strip()
            sample_name = filename.split(".")[0]
            extension = filename[len(sample_name):]
            url= "remote_UCDavis_GoogleDr:projects/SV/results/" + file_path
            destination = os.path.join(output_folder, "results/inputs/")
            command = f"rclone -v --copy-links copy {url} {destination}"
    #        command = f"touch {destination}"
            print(command)
            relative_path = os.path.join("results/inputs/", filename)
            subsample_table[sample].append(relative_path)
            try:
                run_command(command)
                print("Success!")
            except RuntimeError as e:
                print(f"Caught an exception: {e}")
    write_subsample_table(output_folder , subsample_table)

if __name__ == '__main__':
    main()
