"""
Check log.snakefile.txt with ease.
Usage:
python logSnakemakeParser.py log.snakefile.txt
"""

import sys
import re
filename = sys.argv[1]

# Keep track of
lineMetaStart = 0
lineMetaEnd   = 0
lineRule      = []
lineFinish    = []
ruleDict      = {}
# ruleDict structure:
#    {
#   ruleName1: [totalJobNum, [list of submitted id], [list of finished job id]]
#
#}

# Parse out line info for meta info and each rule info
with open(filename) as f:
    for i, line in enumerate(f):
        # if (re.search("^	count	jobs$", line)):
        if line.split() == ["job", "count", "min", "threads", "max", "threads"]:
            lineMetaStart = i + 1
        if (line == "\n"):
            if lineMetaEnd == 0:
                lineMetaEnd   = i
            continue
        if re.search("^rule.*:$", line) or re.search("^localrule.*:$", line):
            if not line.startswith("localrule targets:"):
                lineRule.append(i)
        if re.search("^Finished job .*.$", line):
            lineFinish.append(i)

#  For each rule, parse out total jobs and save in dict:
with open(filename) as f:
    for i, line in enumerate(f):
        if i in range(lineMetaStart + 1, lineMetaEnd - 1):
            tem = line.split()
            jobNum = tem[1]
            ruleName = tem[0]
            # print(line, i, lineMetaStart)
            if not ruleName in ruleDict:
                ruleDict[ruleName] = [jobNum, [], []]

# exit()
# print(ruleDict)

# For each rule, parse out submitted job IDs:
fileList = [line for line in open(filename)]
with open(filename) as f:
    for i, line in enumerate(f):
        if i in lineRule:
            _ruleName = line.split()[1][:-1]
            for j in range(i, i + 6):
                if (re.search('^    jobid: ', fileList[j])):
                    _jobID = fileList[j].split()[1]
                    if not _jobID in ruleDict[_ruleName][1]:
                        ruleDict[_ruleName][1].append(_jobID)
                        break
# For each rule, parse out finished job IDs:
with open(filename) as f:
    for i, line in enumerate(f):
        if i in lineFinish:
            _jobIDFinished = line.split()[2][:-1]
            for key in ruleDict:
                if _jobIDFinished in ruleDict[key][1]:
                    if _jobIDFinished not in ruleDict[key][2]:
                        ruleDict[key][2].append(_jobIDFinished)
                        break

# for key in ruleDict:
#     print(key, ruleDict[key])
for key in ruleDict:
    jt = ruleDict[key][0] # total job number
    js = len(ruleDict[key][1]) # submitted job number
    jf = len(ruleDict[key][2]) # finished job number
    if key == "targets":
        print("Localrule " + (key + ":").ljust(30, " ") + "total -> " + str(jt).ljust(5, " ") + "finished -> " + str(jf).ljust(5, " ") + "percent -> " + str(jf/float(jt)*100).ljust(5, " ") + "%")
    else:
        print("Rule " + (key + ":").ljust(35, " ") + "total -> " + str(jt).ljust(5, " ") + "finished -> " + str(jf).ljust(5, " ") + "percent -> " + str(jf/float(jt)*100).ljust(5, " ") + "%")
