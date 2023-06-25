#<陳羽暉, 111753156>

import argparse
import numpy as np
import math

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="inputpath, pls enter a path")
parser.add_argument("--score", type=str, help="inputscorepath, pls enter a path")
parser.add_argument("--aln", type=str, help="inputaln, pls enter global or local")
parser.add_argument("--gap_open", type=int, help="input gap open score, pls enter a score")
parser.add_argument("--gap_extend", type=int, help="input gap extend score, pls enter a score")
parser.add_argument("--output", type=str, help="outputpath, pls enter a path")
args = parser.parse_args()

file = open(args.score, "r")
data = file.readlines()
file.close()

pamx = np.zeros((24, 24))
for i in range(10, 34):
    s = data[i].rstrip()
    rowData = s.split()
    del rowData[0]
    for j in range(len(rowData)):
        pamx[i - 10][j] = int(rowData[j])
aminoAcid = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']

file = open(args.input, 'r')
ls = []
name = []
for line in file:
    if line.startswith('>'):
        name.append(line)
    if not line.startswith('>'):
        ls.append(line.replace('\n', ''))
file.close()

direction = ['Left', 'Upper Left', 'Up']
gap_open_score = args.gap_open
gap_extend_score = args.gap_extend

dp_score_global = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))
dp_score_global[0][0] = 1
dp_direction_global = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))

ix = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))
ix[0][0] = gap_open_score
iy = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))
iy[0][0] = gap_open_score

for i in range(1, len(ls[0]) + 1):
    dp_score_global[i][0] = -float('inf')
    ix[i][0] = ix[i-1][0] + gap_extend_score
    iy[i][0] = -float('inf')
for j in range(1, len(ls[1]) + 1):
    dp_score_global[0][j] = -float('inf')
    ix[0][j] = -float('inf')
    iy[0][j] = iy[0][j-1] + gap_extend_score

for i in range(1, len(ls[0]) + 1):
    for j in range(1, len(ls[1]) + 1):
        aminoIndex1 = aminoAcid.index(ls[0][i - 1])
        aminoIndex2 = aminoAcid.index(ls[1][j - 1])
        m1 = dp_score_global[i-1][j-1] + pamx[aminoIndex1][aminoIndex2]
        m2 = ix[i-1][j-1] + pamx[aminoIndex1][aminoIndex2]
        m3 = iy[i-1][j-1] + pamx[aminoIndex1][aminoIndex2]
        if (max(m1, m2, m3) == m1):
            dp_score_global[i][j] = m1
        elif (max(m1, m2, m3) == m2):
            dp_score_global[i][j] = m2
        else:
            dp_score_global[i][j] = m3

        ix1 = dp_score_global[i-1][j] + gap_open_score
        ix2 = ix[i-1][j] + gap_extend_score
        if (max(ix1, ix2) == ix1):
           ix[i][j] =  ix1
        else:
            ix[i][j] = ix2

        iy1 = dp_score_global[i][j-1] + gap_open_score
        iy2 = ix[i][j-1] + gap_extend_score
        if (max(iy1, iy2) == iy1):
           iy[i][j] =  iy1
        else:
            iy[i][j] = iy2

for i in range(1, len(ls[0]) + 1):
    for j in range(1, len(ls[1]) + 1):
        if (max(dp_score_global[i][j], ix[i][j], iy[i][j]) == dp_score_global[i][j]):
            dp_direction_global[i][j] = 2
        elif (max(dp_score_global[i][j], ix[i][j], iy[i][j]) == ix[i][j]):
            dp_direction_global[i][j] = 3
        elif (max(dp_score_global[i][j], ix[i][j], iy[i][j]) == iy[i][j]):
            dp_direction_global[i][j] = 1
#print(dp_score_global[33][37], " ", ix[34][37], " ", iy[34][37], " ", dp_direction_global[34][37], " ", ls[0][34], " ", ls[1][37])

dp_score_local = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))
dp_score_local[0][0] = 1
dp_direction_local = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))

ix_local = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))
ix_local[0][0] = gap_open_score
iy_local = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))
iy_local[0][0] = gap_open_score

for i in range(1, len(ls[0]) + 1):
    dp_score_local[i][0] = 0
    ix_local[i][0] = ix_local[i-1][0] + gap_extend_score
    iy_local[i][0] = 0
    if (ix_local[i][0] < 0):
        ix_local[i][0] = 0
for j in range(1, len(ls[1]) + 1):
    dp_score_local[0][j] = 0
    ix_local[0][j] = 0
    iy_local[0][j] = iy[0][j-1] + gap_extend_score
    if (iy_local[i][0] < 0):
        iy_local[i][0] = 0

for i in range(1, len(ls[0]) + 1):
    for j in range(1, len(ls[1]) + 1):
        aminoIndex1 = aminoAcid.index(ls[0][i - 1])
        aminoIndex2 = aminoAcid.index(ls[1][j - 1])
        m1 = dp_score_local[i-1][j-1] + pamx[aminoIndex1][aminoIndex2]
        m2 = ix_local[i-1][j-1] + pamx[aminoIndex1][aminoIndex2]
        m3 = iy_local[i-1][j-1] + pamx[aminoIndex1][aminoIndex2]
        if (max(m1, m2, m3) == m1):
            dp_score_local[i][j] = m1
        elif (max(m1, m2, m3) == m2):
            dp_score_local[i][j] = m2
        else:
            dp_score_local[i][j] = m3
        if (dp_score_local[i][j] < 0):
            dp_score_local[i][j] == 0

        ix1 = dp_score_local[i-1][j] + gap_open_score
        ix2 = ix_local[i-1][j] + gap_extend_score
        if (max(ix1, ix2) == ix1):
           ix_local[i][j] =  ix1
        else:
            ix_local[i][j] = ix2
        if (ix_local[i][j] < 0):
            ix_local[i][j] = 0

        iy1 = dp_score_local[i][j-1] + gap_open_score
        iy2 = ix_local[i][j-1] + gap_extend_score
        if (max(iy1, iy2) == iy1):
           iy_local[i][j] =  iy1
        else:
            iy_local[i][j] = iy2
        if (iy_local[i][j] < 0):
            iy_local[i][j] = 0

for i in range(1, len(ls[0]) + 1):
    for j in range(1, len(ls[1]) + 1):
        if (max(dp_score_local[i][j], ix_local[i][j], iy_local[i][j]) == dp_score_local[i][j]):
            dp_direction_local[i][j] = 2
        elif (max(dp_score_local[i][j], ix_local[i][j], iy_local[i][j]) == ix_local[i][j]):
            dp_direction_local[i][j] = 3
        elif (max(dp_score_local[i][j], ix_local[i][j], iy_local[i][j]) == iy_local[i][j]):
            dp_direction_local[i][j] = 1

if (args.aln == 'global'):
    aminoAcid1 = []
    aminoAcid2 = []
    row = len(ls[0])
    column = len(ls[1])

    traceback = dp_direction_global[row][column]

    while (traceback != 0):
        if traceback == 3:
            aminoAcid1.append(ls[0][row-1])
            aminoAcid2.append('-')
            row -= 1
        elif traceback == 2:
            aminoAcid1.append(ls[0][row-1])
            aminoAcid2.append(ls[1][column-1])
            row -= 1
            column -= 1
            traceback = dp_direction_global[row][column]
        elif traceback == 1:
            aminoAcid1.append('-')
            aminoAcid2.append(ls[1][column-1])
            column -= 1
        traceback = dp_direction_global[row][column]

    aminoAcid1.reverse()
    aminoAcid2.reverse()
    result = [aminoAcid1, aminoAcid2]

    file = open(args.output, 'w')
    lines = []
    for i in range(len(result)):
        lines.append(name[i])
        a = ''.join(result[i]) + "\n"
        lines.append(a)
    file.writelines(lines)
    file.close()

elif (args.aln == 'local'):
    maxScore = 0
    for i in range(len(ls[0]) + 1):
        for j in range(len(ls[1]) + 1):
            if max(dp_score_local[i][j], ix_local[i][j], iy_local[i][j]) > maxScore:
                maxScore = max(dp_score_local[i][j], ix_local[i][j], iy_local[i][j])
    maxIndex = []
    for i in range(len(ls[0]) + 1):
        for j in range(len(ls[1]) + 1):
            if max(dp_score_local[i][j], ix_local[i][j], iy_local[i][j]) == maxScore:
                maxIndex.append((i, j))
    #print(maxIndex)
    result = []
    for i in range(len(maxIndex)):
        aminoAcid1 = []
        aminoAcid2 = []
        (row, column) = maxIndex[i]

        traceback = dp_direction_local[row][column]

        while (traceback != 0):
            if traceback == 3:
                aminoAcid1.append(ls[0][row-1])
                aminoAcid2.append('-')
                row -= 1
            elif traceback == 2:
                aminoAcid1.append(ls[0][row-1])
                aminoAcid2.append(ls[1][column-1])
                row -= 1
                column -= 1
            elif traceback == 1:
                aminoAcid1.append('-')
                aminoAcid2.append(ls[1][column-1])
                column -= 1
            traceback = dp_direction_local[row][column]
        aminoAcid1.reverse()
        aminoAcid2.reverse()
        result.append(aminoAcid1)
        result.append(aminoAcid2)
    
    finalResult = []
    maxLen = 0
    for i in range(len(result)):
        if len(result[i]) > maxLen:
            finalResult = []
            finalResult.append(result[i])
            maxLen = len(result[i])
        elif len(result[i]) == maxLen:
            finalResult.append(result[i])
    
    file = open(args.output, 'w')
    lines = []
    for i in range(0, len(finalResult), 2):
        lines.append(name[0])
        a = ''.join(finalResult[i]) + "\n"
        lines.append(a)
        lines.append(name[1])
        b = ''.join(finalResult[i+1]) + "\n"
        lines.append(b)
    file.writelines(lines)
    file.close()


#for i in range(0, len(ls[0]) + 1):
#    for j in range(0, len(ls[1]) + 1):
#        print(" ", dp_direction_global[i][j], ": ", i, ", ", j, " ", dp_score_global[i][j], " | ", ix[i][j], " | ", iy[i][j], end=' || ')
#    print('\n')

#for i in range(0, len(ls[0]) + 1):
#    for j in range(0, len(ls[1]) + 1):
#        print(" ", dp_direction_global[i][j], end=' ')
#    print('\n')