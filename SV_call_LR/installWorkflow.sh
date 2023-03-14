directory=$1
mkdir -p $directory
cp config.yaml happy.yaml  project_config.yaml  sample_table.csv subsample_table.csv Snakefile $directory
