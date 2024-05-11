import pandas as pd

fname = snakemake.input[0]  # config['META']
# fname = 'meta.csv'
if fname.endswith('csv'):
    df = pd.read_csv(fname)

df.columns = ['sample_name','group','batch']
df['BAM_file'] = 'mapped_reads/' + df['sample_name'] + '.bam'
df['salmon_quant_file'] = 'salmon/' + df['sample_name'] + '/quant.sf'

df.to_csv(snakemake.output[0],  index=False)
