import os

def are_files_non_empty(filenames):
  """Checks if all files in the given list are non-empty.

  Args:
    filenames: A list of file paths.

  Returns:
    True if all files are non-empty, False otherwise.
  """

  for filename in filenames:
    if os.path.getsize(filename) == 0:
      return False
  return True

# Example usage:
file_list = snakemake.input[[0]]
print(file_list)
if are_files_non_empty(file_list):
  print("All files are non-empty")
else:
  print("At least one file is empty")
