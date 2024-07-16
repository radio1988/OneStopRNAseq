import os
import sys


def are_files_non_empty(filenames):
  """Checks if all files in the given list are non-empty.

  Args:
    filenames: A list of file paths.

  Returns:
    True if all files are non-empty, False otherwise.
  """

  for filename in filenames:
    print(filename, os.path.getsize(filename), output)
    if os.path.getsize(filename) < 100:
      return False
  return True

# Example usage:
output = open(snakemake.output[0], "w")

file_list = snakemake.input
print(file_list, file=output)
if are_files_non_empty(file_list):
  print("All files are non-empty", output)
else:
  print("At least one file is empty", output)
  sys.exit("At least trimmed fastq.gz is too small")

output.close()
