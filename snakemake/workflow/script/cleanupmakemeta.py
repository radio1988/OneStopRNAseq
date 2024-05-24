import sys

# Define the log file path (replace with your desired path)
log_file = snakemake.log[0]

# Open the log file for writing in append mode
with open(log_file, "a") as log:
    # Redirect standard error and standard output to the log file
    sys.stderr = log
    sys.stdout = log

    import pandas as pd

    fname = snakemake.input[0]  # config['META']
    if fname.endswith('csv'):
        df = pd.read_csv(fname)

    df.columns = ['sample_name', 'group', 'batch']
    df['BAM_file'] = 'mapped_reads/' + df['sample_name'] + '.bam'
    df['salmon_quant_file'] = 'salmon/A/' + df['sample_name'] + '/quant.sf'

    df.to_csv(snakemake.output[0], index=False)
