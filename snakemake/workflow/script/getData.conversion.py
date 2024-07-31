import re

def extract_two_words_after_pattern(text, pattern):
  """Extracts the two words after the specified pattern in a text string.

  Args:
    text: The input text string.
    pattern: The pattern to search for.

  Returns:
    A tuple containing the two words after the pattern, or None if not found.

  e.g. input:
  bsub -q long -n 1 -W 12:00 -R rusage[mem=500] -o /home/kai.hu-umw/pi/hira.goel-umw/OneStopRNAseq/users/651/GSE190667_CEBPGvsCHOP_02.16.55-07.30.2024/fastq/bsub_log_SRX13375781.txt '( cd /home/kai.hu-umw/pi/hira.goel-umw/OneStopRNAseq/users/651/GSE190667_CEBPGvsCHOP_02.16.55-07.30.2024/fastq && python ../script/getData.py SRX13375781 GSE190667_CHOP3d5 /home/kai.hu-umw/pi/hira.goel-umw/OneStopRNAseq/users/651/GSE190667_CEBPGvsCHOP_02.16.55-07.30.2024/fastq pair ) >> /home/kai.hu-umw/pi/hira.goel-umw/OneStopRNAseq/users/651/GSE190667_CEBPGvsCHOP_02.16.55-07.30.2024/log_ssh2/SRX13375781_getDataPy.txt' &
  """
  match = re.search(rf"{pattern}\s+(\w+)\s+(\w+)\s+(\S+)\s+(\w+)", text)
  if match:
    return match.groups()
  else:
    return None

fname = 'getData.sh'
pattern = '../script/getData.py'

print("""
# bsub -W 144:00 'bash download.sh'
set -e
source /pi/mccb-umw/shared/conda/miniconda3/etc/profile.d/conda.sh
conda activate /pi/mccb-umw/shared/OneStopRNAseq/conda/osr-base
module load sratoolkit/3.0.0
module load perl/5.36.0_AS
""")

with open(fname, 'r') as file:
    for line in file:
        if pattern in line:
            accession, name, path, type = extract_two_words_after_pattern(line, pattern)
            if accession is not None:
                print('python getData.py {} {} . {}'.format(accession, name, type))

# example output:
# bsub -W 144:00 'bash download.sh'
# set -e
# source /pi/mccb-umw/shared/conda/miniconda3/etc/profile.d/conda.sh
# conda activate /pi/mccb-umw/shared/OneStopRNAseq/conda/osr-base
# module load sratoolkit/3.0.0
# module load perl/5.36.0_AS
# python getData.py SRX13375766 GSE190667_ATF33d5 . pair
# python getData.py SRX13375767 GSE190667_ATF33d6 . pair
# ...
