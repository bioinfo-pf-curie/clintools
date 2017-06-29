#! /usr/bin/env python
# -*- coding: utf8 -*-

# This file is part of table_maker.py

# Copyright Institut Curie 2014

# This software is a computer program whose purpose is to MaxEntScan scores.

# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

import optparse, sys, re, string, subprocess, tempfile, os
import clinTools

parser = optparse.OptionParser(usage="%prog [options] conf_file ref_genome_file",
    description="Generate a report of number of reads supporting each base at given regions. This script work on mpileup and bam files depending on files extensions in conf file (.pileup or .bam). If you use bam files, it is highly recommended to use a bed file with the -b option. The configuration must contain \"sample_name\\tpath_to_file(.pileup or .bam)\" in each line. Output positions are 1-based")
parser.add_option("-r", "--ratio", default=False, action="store_true",
    help="If set allelic ratio are reported in output"
)
parser.add_option("-n", "--names", default=False, action="store_true",
    help="Print names from bed in output"
)
parser.add_option("-b", "--bed_file", default=None, type=str,
    help="Bed file to use for intersection (0-based positions)"
)
parser.add_option("-q", "--mapQ", default=0,
    help="-q option for mpileup. mapQ threshold (default : 0). This option is useful only if you use bam files as input."
)
parser.add_option("-Q", "--BAQ", default=0,
    help="-Q option for mpileup. base quality threshold (default : 0). This option is useful only if you use bam files as input."
)
parser.add_option("-s", "--keep-strand", default=False, action="store_true",
    help="If set base counts for both strands are separated"
)
(options, args) = parser.parse_args()

if len(args) != 2:
    parser.print_help()
    exit(1)

clinTools.checkSamtoolsVersion('1.1')

confFile = open(args[0], 'r')
reference = args[1]

bedFileName = options.bed_file

#stderrFile = open('/dev/null', 'a') 
stderrFile = sys.stderr 

if options.keep_strand:
    if options.ratio:
        headerLineList = ["barcode", "chromosome", "position", 'reference', 'depth', 'A+', 'A+_ratio', 'A-', 'A-_ratio', 'T+', 'T+_ratio', 'T-', 'T-_ratio', 'C+', 'C+_ratio', 'C-', 'C-_ratio', 'G+', 'G+_ratio', 'G-', 'G-_ratio', 'N', 'N_ratio', 'Ins', 'Del']
    else:
        headerLineList = ["barcode", "chromosome", "position", 'reference', 'depth', 'A+', 'A-', 'T+', 'T-', 'C+', 'C-', 'G+', 'G-', 'N', 'Ins', 'Del']
else:
    if options.ratio:
        headerLineList = ["barcode", "chromosome", "position", 'reference', 'depth', 'A', 'A_ratio', 'T', 'T_ratio', 'C', 'C_ratio', 'G', 'G_ratio', 'N', 'N_ratio', 'Ins', 'Del']
    else:
        headerLineList = ["barcode", "chromosome", "position", 'reference', 'depth', 'A', 'T', 'C', 'G', 'N', 'Ins', 'Del']

if bedFileName != None:
    if options.names:
        headerLineList.append('Name')
        bedFile = open(bedFileName, 'r')
        bedDict = dict()
        chromList = list()
        for line in bedFile:
            lineList = line.rstrip().split('\t')
            if lineList[0] in bedDict.keys():
                bedDict[lineList[0]].append({'start' : str(int(lineList[1]) + 1), 'end' : str(int(lineList[2]) + 1), 'name' : lineList[3]})
            else:
                bedDict[lineList[0]] = [{'start' : str(int(lineList[1]) + 1), 'end' : str(int(lineList[2]) + 1), 'name' : lineList[3]}]
                chromList.append(lineList[0])
    else:
        bedFile = open(bedFileName, 'r')
        bedDict = dict()
        chromList = list()
        for line in bedFile:
            lineList = line.rstrip().split('\t')
            # convert start position to 1-based but not end to have the same behaviour as samtools mpileup -l
            if lineList[0] in bedDict.keys():
                bedDict[lineList[0]].append({'start' : str(int(lineList[1]) + 1), 'end' : str(int(lineList[2]))})
            else:
                bedDict[lineList[0]] = [{'start' : str(int(lineList[1]) + 1), 'end' : str(int(lineList[2]))}]
                chromList.append(lineList[0])


# print header
print '\t'.join(headerLineList)

tempFileName = None

for confLine in confFile:
    try:
        barcode, fileName = confLine.rstrip().split('\t')

        if re.search('.bam$', fileName):
            clinTools.checkSamtoolsVersion('1.1')
            tempPileupFile = tempfile.NamedTemporaryFile(mode="a+b", delete = False)
            tempFileName = tempPileupFile.name
            cmdList = ["samtools", "mpileup", "-AB", "-f", reference, "-q", str(options.mapQ), "-Q", str(options.BAQ), "-d 1000000"]
            if bedFileName != None:
                for chrom in chromList:
                    for region in bedDict[chrom]:
                        cmdList1 = list(cmdList)

                        cmdList1.append("-r " + chrom + ":" + region['start'] + "-" + region['end'])
                        cmdList1.append(fileName)

                        cmd = ' '.join(cmdList1)
                        sp = subprocess.Popen(cmd, stdout = tempPileupFile, stderr=stderrFile, shell=True)
                        sp.wait()
            else:
                cmdList.append(fileName)
                cmd = ' '.join(cmdList)
                sp = subprocess.Popen(cmd, stdout = tempPileupFile, stderr=stderrFile, shell=True)
                sp.wait()

            tempPileupFile.close()
            pileupFile = open(tempPileupFile.name, 'r')
        elif re.search('.pileup', fileName):
            if bedFileName != None:
                tempPileupFile = tempfile.NamedTemporaryFile(mode = 'w', delete = False)
                tempFileName = tempPileupFile.name
                # intersect bed with pileup
                pileupFile = open(fileName, 'r')
                for linePileup in pileupFile:
                    chrom, pos = linePileup.rstrip().split('\t')[:2]
                    if chrom in chromList:
                        for region in bedDict[chrom]:
                            if int(region['start']) <= int(pos) <= int(region['end']):
                                tempPileupFile.write(linePileup)
                                break

                tempPileupFile.close()
                pileupFile = open(tempFileName, 'r')
            else:
                pileupFile = open(fileName, 'r')
        else:
            raise Exception("invalid file extension. usable extension are .pileup and .bam")

        # create storing variable for dels --> useful to report deletions at the right position
        delStore = None

        for line in pileupFile:
            pileupDict = clinTools.parsePileupLine(line)
            pos = pileupDict['pos']

            outLineList = [barcode, pileupDict['chrom'], pos, pileupDict['ref'], str(pileupDict['depth'])]

            if options.keep_strand:
                if options.ratio:
                    if pileupDict['depth'] != 0:
                        outLineList.append(pileupDict['A'][0])
                        outLineList.append("%.2f" % ((float(pileupDict['A'][0])/float(pileupDict['depth']))*100.0))
                        outLineList.append(pileupDict['A'][1])
                        outLineList.append("%.2f" % ((float(pileupDict['A'][1])/float(pileupDict['depth']))*100.0))
                        outLineList.append(pileupDict['T'][0])
                        outLineList.append("%.2f" % ((float(pileupDict['T'][0])/float(pileupDict['depth']))*100.0))
                        outLineList.append(pileupDict['T'][1])
                        outLineList.append("%.2f" % ((float(pileupDict['T'][1])/float(pileupDict['depth']))*100.0))
                        outLineList.append(pileupDict['C'][0])
                        outLineList.append("%.2f" % ((float(pileupDict['C'][0])/float(pileupDict['depth']))*100.0))
                        outLineList.append(pileupDict['C'][1])
                        outLineList.append("%.2f" % ((float(pileupDict['C'][1])/float(pileupDict['depth']))*100.0))
                        outLineList.append(pileupDict['G'][0])
                        outLineList.append("%.2f" % ((float(pileupDict['G'][0])/float(pileupDict['depth']))*100.0))
                        outLineList.append(pileupDict['G'][1])
                        outLineList.append("%.2f" % ((float(pileupDict['G'][1])/float(pileupDict['depth']))*100.0))
                        outLineList.append(pileupDict['countN'])
                        outLineList.append("%.2f" % ((float(pileupDict['countN'])/float(pileupDict['depth']))*100.0))
                    else:
                        outLineList.extend(pileupDict['A'])
                        outLineList.append("%.2f" % 0)
                        outLineList.extend(pileupDict['T'])
                        outLineList.append("%.2f" % 0)
                        outLineList.extend(pileupDict['C'])
                        outLineList.append("%.2f" % 0)
                        outLineList.extend(pileupDict['G'])
                        outLineList.append("%.2f" % 0)
                        outLineList.append(pileupDict['countN'])
                        outLineList.append("%.2f" % 0)
                else:
                    outLineList.extend(pileupDict['A'])
                    outLineList.extend(pileupDict['T'])
                    outLineList.extend(pileupDict['C'])
                    outLineList.extend(pileupDict['G'])
                    outLineList.append(pileupDict['countN'])

                # add ins string
                insString = None
                for key, value in zip(pileupDict['ins'].keys(), pileupDict['ins'].values()):
                    if insString == None:
                        insString = "{seq}:{countFw},{countRv}".format(seq = key, countFw = value[0], countRv = value[1])
                    else:
                        insString += ";{seq}:{countFw},{countRv}".format(seq = key, countFw = value[0], countRv = value[1])
                outLineList.append(insString)

                # add del string
                delString = None
                storeString = None
                for key, value in zip(pileupDict['del'].keys(), pileupDict['del'].values()):
                    if key == '*':
                        addString = "{seq}:{count}".format(seq = key, count = sum(value))
                        if delString == None:
                            delString = addString
                        else:
                            delString += ";" + addString
                    else:
                        addString = "{seq}:{countFw},{countRv}".format(seq = key, countFw = value[0], countRv = value[1])
                        if storeString == None:
                            storeString = addString 
                        else:
                            storeString += ";" + addString

                if delStore != None:
                    if delString == None:
                        delString = delStore
                    else:
                        delString += ';' + delStore
                outLineList.append(delString)
                delStore = storeString

            else:
                if options.ratio:
                    if pileupDict['depth'] != 0:
                        outLineList.append(sum(pileupDict['A']))
                        outLineList.append("%.2f" % ((float(sum(pileupDict['A']))/float(pileupDict['depth']))*100.0))
                        outLineList.append(sum(pileupDict['T']))
                        outLineList.append("%.2f" % ((float(sum(pileupDict['T']))/float(pileupDict['depth']))*100.0))
                        outLineList.append(sum(pileupDict['C']))
                        outLineList.append("%.2f" % ((float(sum(pileupDict['C']))/float(pileupDict['depth']))*100.0))
                        outLineList.append(sum(pileupDict['G']))
                        outLineList.append("%.2f" % ((float(sum(pileupDict['G']))/float(pileupDict['depth']))*100.0))
                        outLineList.append(pileupDict['countN'])
                        outLineList.append("%.2f" % ((float(pileupDict['countN'])/float(pileupDict['depth']))*100.0))
                    else:
                        outLineList.append(sum(pileupDict['A']))
                        outLineList.append("%.2f" % 0)
                        outLineList.append(sum(pileupDict['T']))
                        outLineList.append("%.2f" % 0)
                        outLineList.append(sum(pileupDict['C']))
                        outLineList.append("%.2f" % 0)
                        outLineList.append(sum(pileupDict['G']))
                        outLineList.append("%.2f" % 0)
                        outLineList.append(pileupDict['countN'])
                        outLineList.append("%.2f" % 0)
                else:
                    outLineList.append(sum(pileupDict['A']))
                    outLineList.append(sum(pileupDict['T']))
                    outLineList.append(sum(pileupDict['C']))
                    outLineList.append(sum(pileupDict['G']))
                    outLineList.append(pileupDict['countN'])

                # add ins string
                insString = 0
                for key, value in zip(pileupDict['ins'].keys(), pileupDict['ins'].values()):
                    if insString == 0:
                        insString = "{seq}:{count}".format(seq = key, count = sum(value))
                    else:
                        insString += ";{seq}:{count}".format(seq = key, count = sum(value))
                outLineList.append(insString)

                # add del string
                delString = "0"
                storeString = "0"
                for key, value in zip(pileupDict['del'].keys(), pileupDict['del'].values()):
                    addString = "{seq}:{count}".format(seq = key, count = sum(value))
                    if key == '*':
                        if delString == "0":
                            delString = addString
                        else:
                            delString += ";" + addString
                    else:
                        if storeString == "0":
                            storeString = addString 
                        else:
                            storeString += ";" + addString

                if delStore != None:
                    if delString == "0":
                        delString = delStore
                    else:
                        delString += ';' + delStore
                outLineList.append(delString)
                delStore = storeString

            if bedFileName != None:
                if options.names:
                    namesList = list()
                    for elt in bedDict[pileupDict['chrom']]:
                        if int(elt['start']) <= int(pos) <= int(elt['end']):
                            namesList.append(elt['name'])

                    outLineList.append(','.join(namesList))

            print '\t'.join([str(elt) for elt in outLineList])
    finally:
        if tempFileName:
            os.remove(tempFileName)
