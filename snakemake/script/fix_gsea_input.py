import sys

print("function: make all lower case letters in txt file upper case, then save back to the same file")
print("usage: capslock.py file.txt")

fname = sys.argv[1]

inputFile = open(fname, 'r')
content = inputFile.read()

# search for lower letters
contain_lower = False
for line in content:
    if any(letter.islower() for letter in line):
        contain_lower = True
        print("log: lower case letters found, updating", fname)
        content=content.upper()
        break

# search for spaces
contain_space = False
for line in content:
    if ' ' in line:
        contain_space = True
        print("log: spaces found in lines, updating", fname)
        content=content.replace(" ", "") # remove space
        break

# convert to upper 
if contain_lower or contain_space:
    with open(fname, 'w') as outputFile:
        outputFile.write(content)
else:
    print("lower case nor spaces found in file, not updating", fname)

