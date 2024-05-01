import re
import numpy as np

alist_filename = 'PG(3,11).alist'
input_filename = '../Outputs/OutputPG(3,11).txt'


def parse_GAP_header(header):
    # Dirty trick of [\S\s] is space or not space
    p = re.compile('(\d+)b x (\d+)[\S\s]+b is (\d+)')
    arrStrings = p.findall(header)[0]
    return int(arrStrings[0]), int(arrStrings[1]), int(arrStrings[2])


def splitData(data):
    header = ""
    startHrep = False
    dataH = ""
    dataP = ""
    for singleLine in data:
        if "Pbin is" in singleLine:
            break
        elif startHrep:
            dataH += singleLine.strip()
        elif "Hrep is" in singleLine:
            startHrep = True
        else:
            header += singleLine.strip()
    return header, dataH, dataP


def read_GAP_output(filename):
    myFile = open(filename, 'r')
    data = myFile.readlines()
    myFile.close()
    header, data, dataP = splitData(data)
    numRows, numCols, b = parse_GAP_header(header)
    return numRows, numCols, b, readHrep(data, numRows, numCols, b)


def readHrep(data, numRows, numCols, b):
    p = re.compile('\[(?:[^\]\[]+)*\]')
    arrStrings = p.findall(data)
    Hrep = [[None] * numCols] * numRows
    for i in np.arange(numRows):
        Hrep[i] = [None] * numCols  # so we make a different array
        for j in np.arange(numCols):
            singleCell = arrStrings[i * numCols + j]
            # take out parenthesis and spaces
            cleanString = singleCell[1: -1].replace(' ', '')
            arrNums = cleanString.split(',')  # we split by commas
            if arrNums[0] == '':
                Hrep[i][j] = []
            else:
                Hrep[i][j] = np.array([int(s) for s in arrNums])
    return Hrep


def get_nonzeros_row(row, sft, b):
    cter = 0
    nonzeros = []
    for block in row:
        if len(block) > 0:
            sftedBlock = (block + sft) % b
            for el in sftedBlock:
                nonzeros.append(cter + el)
        cter += b
    return nonzeros


def get_nonzeros_col(col, sft, b):
    cter = 0
    nonzeros = []
    for block in col:
        if len(block) > 0:
            sftedBlock = (b - block + sft) % b
            for el in sftedBlock:
                nonzeros.append(cter + el)
        cter += b
    return nonzeros


def write_alist_file(filename, Hrep, numRows, numCols, b):
    """
    This function writes an alist file for the parity check
    matrix. The format of alist files is described at:
    http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html
    """
    with open(filename, 'w') as myfile:
        numRowsb = numRows * b
        numColsb = numCols * b
        tempstring = repr(numColsb) + ' ' + repr(numRowsb) + '\n'
        myfile.write(tempstring)
        tempstring1 = ''
        tempstring2 = ''
        maxRowWeight = 0
        for rowNum in np.arange(numRows):
            rowHrep = Hrep[rowNum]
            for sft in np.arange(b):
                rowNumb = rowNum * b + sft
                nonzeros = get_nonzeros_row(rowHrep, sft, b)
                rowWeight = len(nonzeros)
                maxRowWeight = max(rowWeight, maxRowWeight)
                tempstring1 += repr(rowWeight) + ' '
                for index in nonzeros:
                    tempstring2 += repr(index + 1) + ' '
                tempstring2 += '\n'
            tempstring1 += '\n'
        tempstring3 = ''
        tempstring4 = ''
        maxColWeight = 0
        for colNum in np.arange(numCols):
            colHrep = []
            for idx in np.arange(numRows):
                colHrep.append(Hrep[idx][colNum])
            for sft in np.arange(b):
                colNumb = colNum * b + sft
                nonzeros = get_nonzeros_col(colHrep, sft, b)
                colWeight = len(nonzeros)
                maxColWeight = max(colWeight, maxColWeight)
                tempstring3 += repr(colWeight) + ' '
                for index in nonzeros:
                    tempstring4 += repr(index + 1) + ' '
                tempstring4 += '\n'
            tempstring3 += '\n'
        tempstring = repr(maxColWeight) + ' ' + repr(maxRowWeight) + '\n'
        # write out max column and row weights
        myfile.write(tempstring)
        # write out all of the column weights
        myfile.write(tempstring3)
        # write out all of the row weights
        myfile.write(tempstring1)
        # write out the nonzero indices for each column
        myfile.write(tempstring4)
        # write out the nonzero indices for each row
        myfile.write(tempstring2)


numRows, numCols, b, Hrep = read_GAP_output(input_filename)
write_alist_file(alist_filename, Hrep, numRows, numCols, b)
