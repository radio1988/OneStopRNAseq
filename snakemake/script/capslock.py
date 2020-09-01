import sys

print("function: make all lower case letters in txt file upper case, then save back to the same file")
print("usage: capslock.py file.txt")

fname = sys.argv[1]

inputFile = open(fname, 'r')
content = inputFile.read()
contain_lower = False

for line in content:
    if any(letter.islower() for letter in line):
        contain_lower = True
        break

if contain_lower:
    print("log: lower case letters found, updating", fname)
    with open(fname, 'w') as outputFile:
        outputFile.write(content.upper())
else:
    print("lower case not found, not updating", fname)

