#!/usr/bin/env python3
'''
File: sam_abundance_modified.py
Author: Shixu He
Email: heshixu@genomics.cn
Date: 2018-08-21
----
Calculator for abundance in bowtie/bwa pipeline (cOMG).
----
Version: 0.2
Add output of read-marker map file

Version: 0.3
+Repair the discrepancy between sam abundance calculator and soap abundance calculator.
+Add the funcion about filtering of conserced read list (read aligned to two or more references).
+Change the output format to :
marker_ID,reads_fragment_number,relative_abundance,reads_list
Note that the reads number will be halved for pair end mode (reads_fragment_number = reads_number/2).
+Add the function about filtering of conserved reads (reads aligned to two or more reference).
+Add the function for arguments control.
'''

import sys
import time
import re
import os
import argparse
from collections import OrderedDict
import pysam


class SamRecord(object):
    '''record the key information of each line in sam file'''

    def __init__(self, flag, marker, length, pos, name):
        super(SamRecord, self).__init__()
        try:
            self.readFlag = int(flag)     # the second column of sam file
            # the length of 10th SEQ column of sam file
            self.readLength = int(length)
            self.markerName = marker      # the reference gene name, the 3rd column
            # the location of first aligned base, the 4th column
            self.alignPos = int(pos)
            self.readName = name
        except Exception as e:
            raise e

    def __str__(self):
        if self.isFirstRead():
            return self.readName + "/1"
        else:
            return self.readName + "/2"

    def getReadLength(self):
        return self.readLength

    def getMarkerName(self):
        return self.markerName

    def getFirstLoc(self):
        return self.alignPos

    def isReverse(self):
        if(self.readFlag & 0x10):
            return True
        return False

    def isFirstRead(self):
        if(self.readFlag & 0x40):
            return True
        return False

    def isProperlyAlignPair(self):
        if(self.readFlag & 0x2):
            return True
        return False

    def isPairMode(self):
        if(self.readFlag & 0x1):
            return True
        return False

    def isPrimary(self):
        if(self.readFlag & 0x100):
            return False
        return True

    def getReadName(self):
        if self.isFirstRead():
            return self.readName + "/1"
        else:
            return self.readName + "/2"

    def getReadPrefix(self):
        return self.readName


class RecordMatrix(object):
    '''record for pair read or single read in dict format'''
    '''{readName:
            [
            [read1_SamRecord1, read1_SamRecord2, ...],
            [read2_SamRecord1, Read2_SamRecord2, ...]
            ]
        }

        Update: new format:
        {
            "count1": {readName: <int>}
            "record1": {readName: [samRecord]}
            "count2": {readName: <int>}
            "record2": {readname: []}
        }
    '''

    def __init__(self):
        super(RecordMatrix, self).__init__()
        self.__matrix = {
            "count1": {},
            "record1": {},
            "count2": {},
            "record2": {}
        }

    def addSamRecord(self, readName, samRecord):
        if(samRecord.isFirstRead()):
            self.__matrix["count1"][readName] = self.__matrix["count1"].get(readName, 0) + 1
            if self.__matrix["count1"][readName] > 1:
                return
            self.__matrix["record1"].setdefault(readName, []).append(samRecord)
        else:
            self.__matrix["count2"][readName] = self.__matrix["count2"].get(readName, 0) + 1
            if self.__matrix["count2"][readName] > 1:
                return
            self.__matrix["record2"].setdefault(readName, []).append(samRecord)

    def getRecordCount(self, readName, strand):
        # if(samRecord.isFirstRead()):
        if (strand == 1):
            return int(self.__matrix["count1"].get(readName, 0))
        elif (strand == 2):
            return int(self.__matrix["count2"].get(readName, 0))
        return 0

    def getReadNameList(self):
        return list(set(self.__matrix["record1"].keys()).union(set(self.__matrix["record2"].keys())))

    def getRecord(self, readName, strand):
        if (strand == 1):
            return self.__matrix["record1"].get(readName, [None])[0]
        elif (strand == 2):
            return self.__matrix["record2"].get(readName, [None])[0]


class Sample(object):
    """Sequence Sample, a Sample may include unique/pair or more reads file and SAM result"""

    def __init__(self, sampleName, samFile, insert_size, sam_format, markerMatrixRef, conservedReadInstance, identity=0.95):
        super(Sample, self).__init__()
        self.__sampleName = sampleName
        self.__identity = identity
        self.conservedRead = conservedReadInstance
        self.markerMatrix = markerMatrixRef
        self.sam_format = sam_format

        #self.__samList = self.__ReadSamListFile(samFileList)
        self.__samList = [(samFile, insert_size)]

        # dict: {"frag_num": {marker_name}
        self.__abunMatrix = {
            "frag_num": {},
            "abundance": {},
            "samrec_list": {}
        }

        self.readsNum = {}
        self.abunOfMarker = {}

        self.__totalAbundance = 0
        self.__totalReads = 0

    def __ReadSamListFile(self, samFileList):
        '''SAM-file-list file content format:
            ^filename   insertsize
        '''
        samList = []
        try:
            with open(samFileList, "r") as smlist:
                for line in smlist:
                    tmp = line.split()
                    samList.append((tmp[0], int(tmp[1])))
        except Exception as e:
            raise e
        return samList

    def __readMarkerMatrix(self, refLength):
        '''
            marker gene length file content format:
            ^No.  reference_gene  length    species_name    genus_name  phylum_name
        '''
        #dict: {marker_name: [markerID, gene_length, species_name]}
        pass

    def ReadsCount(self):
        recordMatrix = RecordMatrix()
        for sam in self.__samList:
            # t_rec1=time.time()
            self.__SamFileRead(sam[0], recordMatrix)
            # t_rec2=time.time()-t_rec1
            #print("sam file read time:\t{0}s".format(t_rec2))
            # you use the print here and stderr.write at main, so the log output is strange
        # t_rec1=time.time()
        self.__MarkersCount(sam[1], recordMatrix)
        #print("markers count time:\t{0}s".format(t_rec2))

    def __SamFileRead(self, samfile, recordMatrix):
        if self.sam_format == 'SAM':
            mode = 'r'
        else:
            mode = 'rb'
        sf = pysam.AlignmentFile(samfile, mode)
        for r in sf:
            if not r.is_unmapped:
                tmpNM = r.get_tag("NM:i")
                # infer_query_length
                tmpLen = sum([int(i) for i in re.findall(
                    r"(\d+)(?:M|I|D)", r.cigarstring)])
                if ((1 - tmpNM / tmpLen) >= self.__identity):
                    _readName = re.sub(r"\/[12]$", "", r.query_name)
                    if self.conservedRead.isConserved(_readName):
                        continue
                    _record = SamRecord(
                        r.flag, r.reference_name, r.infer_read_length(), r.reference_start, _readName)
                    recordMatrix.addSamRecord(_readName, _record)

    def __SamLineFilter(self, line):
        '''
        return useful fields of sam record line if satisfy conditions, else return None
        case 1: paired mapped
        case 2: identity
        return : tmp, a list
        '''

        if(re.match(r"^\@", line)):
            return None
        tmp = re.split(r"\s+", line, 11)

        # case 1
        if(int(tmp[1]) & 0x4):
            return None
        # case 2
        try:
            # 'NM:i:' ( "Edit distance to the reference" )
            tmpNM = int(re.search(r"NM:i:(\d+)", tmp[11]).group(1))
            tmpLen = sum([int(i) for i in re.findall(r"(\d+)(?:M|I|D)", tmp[5])])
        except Exception as e:
            raise e
        if((1-tmpNM/tmpLen) < self.__identity):
            return None
        return tmp

    def __MarkersCount(self, ins, recordMatrix):
        assert type(ins) is int, "insert size type error\n"
        for read in recordMatrix.getReadNameList():
            recLen1 = recordMatrix.getRecordCount(read, 1)
            recLen2 = recordMatrix.getRecordCount(read, 2)
            if(recLen1 == 1):
                if(recLen2 == 1):
                    rec1 = recordMatrix.getRecord(read, 1)
                    rec2 = recordMatrix.getRecord(read, 2)
                    if(rec1.getMarkerName() == rec2.getMarkerName()):
                        self.__PairAlignAdd(rec1, rec2, ins)
                    else:
                        self.__SingleAlignAdd(rec1, ins)
                        self.__SingleAlignAdd(rec2, ins)
                else:
                    self.__SingleAlignAdd(recordMatrix.getRecord(read, 1), ins)
            elif(recLen2 == 1):
                self.__SingleAlignAdd(recordMatrix.getRecord(read, 2), ins)
            else:
                continue

    def __SingleAlignAdd(self, rec, ins, discordant=False):
        assert type(ins) is int, "insertsize type error"
        _markerName = rec.getMarkerName()
        if discordant:
            #self.__abunMatrix["samrec_list"].setdefault(_markerName, []).append(rec)
            return True
        if(rec.isReverse()):
            if((self.markerMatrix.getMarkerLength(_markerName)-rec.getFirstLoc()) < (ins+100)):
                self.__abunMatrix["frag_num"][_markerName] = self.__abunMatrix["frag_num"].get(_markerName, 0) + 1
                #self.__abunMatrix["samrec_list"].setdefault(_markerName, []).append(rec)
                return True
        else:
            if(rec.getFirstLoc() < (ins - rec.getReadLength() + 100)):
                self.__abunMatrix["frag_num"][_markerName] = self.__abunMatrix["frag_num"].get(_markerName, 0) + 1
                #self.__abunMatrix["samrec_list"].setdefault(_markerName, []).append(rec)
                return True
        return False

    def __PairAlignAdd(self, rec1, rec2, ins):
        assert type(ins) is int, "insert size type error\n"
        _markerName = rec1.getMarkerName()
        if (rec1.isProperlyAlignPair()):
            self.__abunMatrix["frag_num"][_markerName] = self.__abunMatrix["frag_num"].get(_markerName, 0) + 1
            #self.__abunMatrix["samrec_list"].setdefault(_markerName, []).append(rec1)
            #self.__abunMatrix["samrec_list"].setdefault(_markerName, []).append(rec2)
        else:
            discordant = False
            discordant = self.__SingleAlignAdd(rec1, ins, discordant)
            self.__SingleAlignAdd(rec2, ins, discordant)

    def AbunCalculate(self):
        for _marker in self.markerMatrix.getMarkerIterator():
            self.__totalAbundance += self.__abunMatrix["frag_num"].get(_marker, 0) / self.markerMatrix.getMarkerLength(_marker)
            self.__totalReads += self.__abunMatrix["frag_num"].get(_marker, 0)
            self.__abunMatrix["abundance"][_marker] = self.__abunMatrix["frag_num"].get(_marker, 0) / self.markerMatrix.getMarkerLength(_marker)

    def AbundanceOutput(self, prefix, outputType, suffix=""):
        if suffix:
            suffix = "." + suffix
        try:
            # evaluation
            # marker    reads(pair)_number  relative_abundance  reads_list
            if (outputType == "evaluation") or (outputType == "all"):
                with open(os.path.join(prefix, self.__sampleName + ".evaluation.abundance" + suffix), "w") as abunOut:
                    abunOut.write("marker,reads(pair) number,relative abundance,reads list,\n")
                    for _marker in self.__abunMatrix["frag_num"].keys():
                        abunOut.write("{marker},{frag_num},{rel_abun},{reads_list},\n".format(
                            marker     = _marker,
                            frag_num   = self.__abunMatrix["frag_num"][_marker],
                            rel_abun   = self.__abunMatrix["abundance"][_marker] / self.__totalAbundance,
                            reads_list = "|".join([str(rec) for rec in self.__abunMatrix["samrec_list"][_marker]])))

            # standard
            # ID    reads_pairs reads_proportion    gene_abundance
            elif (outputType == 'standard') or (outputType == "all"):
                with open(os.path.join(prefix, self.__sampleName + ".abundance"  + suffix), "w") as abunOut:
                    abunOut.write("ID\treads_pairs\treads_proportion\tgene_abundance\n")
                    for _marker in self.markerMatrix.getMarkerIterator():
                        abunOut.write("{markerID}\t{reads_num}\t{reads_prop}\t{rel_abun}\n".format(
                            markerID   = self.markerMatrix.getMarkerID(_marker),
                            reads_num  = self.__abunMatrix["frag_num"].get(_marker, 0),
                            reads_prop = self.__abunMatrix["frag_num"].get(_marker, 0)  / self.__totalReads,
                            rel_abun   = self.__abunMatrix["abundance"].get(_marker, 0) / self.__totalAbundance))
            elif (outputType == "non-zero"):
                with open(os.path.join(prefix, self.__sampleName + ".non_zero.abundance" + suffix), "w") as abunOut:
                    abunOut.write(
                        "ID\treads_pairs\treads_proportion\tgene_abundance\n")
                    for _marker in self.markerMatrix.getMarkerIterator():
                        if self.__abunMatrix["frag_num"].get(_marker, 0) == 0:
                            continue
                        abunOut.write("{markerID}\t{reads_num}\t{reads_prop}\t{rel_abun}\n".format(
                            markerID   = self.markerMatrix.getMarkerID(_marker),
                            redes_num  = self.__abunMatrix["frag_num"].get(_marker, 0),
                            reads_prop = self.__abunMatrix["frag_num"].get(_marker, 0)  / self.__totalReads,
                            rel_abun   = self.__abunMatrix["abundance"].get(_marker, 0) / self.__totalAbundance))
            with open(os.path.join(prefix, self.__sampleName + ".abundance.size" + suffix), "w") as abunSizeOut:
                abunSizeOut.write("TotalReads:\t%d\n" % self.__totalReads)

        except IOError as ioe:
            raise ioe
        except Exception as e:
            raise e


class MarkerMatrix(object):
    '''
    Read IGC 9.9M udpate matrix reference file.
    File format:
        markerID    marker_name marker_gene_length  s__species_name g__genus_name   p__phylum_name

    Useful column:
        markerID, marker_name, marker_gene_length
    '''

    def __init__(self, markerMatrixFile):
        super(MarkerMatrix, self).__init__()
        # OrderedDict: {marker_name: {"markerID": <str>, "gene_length": <int>}}
        self.__markerMatrix = OrderedDict()
        self.__scanMarkerMatrixFile(markerMatrixFile)

    def __scanMarkerMatrixFile(self, markerMatrixFile):
        assert os.path.exists(
            markerMatrixFile), "file {0} doesn't exists.\n".format(markerMatrixFile)
        try:
            # format:
            # mgs_id mgs_id_old contig_id contig_id_old contig_name contig_length species genus phylum
            with open(markerMatrixFile, 'r', encoding='utf-8') as _mmf:
                next(_mmf)
                while True:
                    lines = _mmf.readlines(65535)
                    if not lines:
                        break
                    for line in lines:
                        tmp = line.rstrip().split()
                        self.__markerMatrix[tmp[4]] = {
                            "markerID": tmp[2],
                            "gene_length": tmp[5]
                        }
        except IOError as _ioe:
            raise _ioe
        except Exception as _e:
            raise _e

    def getMarkerID(self, markerName):
        return self.__markerMatrix[markerName]["markerID"]

    def getMarkerLength(self, markerName):
        return int(self.__markerMatrix[markerName]["gene_length"])

    def getMarkerIterator(self):
        return self.__markerMatrix.keys()


class ConservedRead(object):
    '''
    Record the ID of conerved read which is proved to have two or more reference during alignment.
    '''

    def __init__(self, conservedReadListFile=None):
        super(ConservedRead, self).__init__()
        self.__conservedReadDict = {}
        self.__readConservedReadList(conservedReadListFile)

    def __readConservedReadList(self, conservedReadListFile):
        if (conservedReadListFile == None) or (conservedReadListFile == ""):
            return
        assert os.path.exists(conservedReadListFile), "file {0} doesn't exists.\n".format(
            conservedReadListFile)
        try:
            with open(conservedReadListFile, 'r', encoding='utf-8') as _crlf:
                while True:
                    lines = _crlf.readlines(65535)
                    if not lines:
                        break
                    for line in lines:
                        self.__conservedReadDict[line.rstrip()] = None
        except IOError as _ioe:
            raise _ioe
        except Exception as _e:
            raise _e
        return

    def isConserved(self, read):
        if re.sub(r"\/[12]$", "", read) in self.__conservedReadDict:
            return True
        return False


def checkArgv():
    '''Parse the command line parameters.'''
    _parser = argparse.ArgumentParser(
        description="The script is used for abundance calculation from alignment of bowtie/bwa.")

    _parser.add_argument(
        "--sample-name",
        required=True,
        help="The sample identity"
    )
    _parser.add_argument(
        "--sam-file",
        required=True,
        help="The paths of SAM file from alignment software"
    )
    _parser.add_argument(
        "--insert-size",
        required=True,
        type=int,
        default=350,
        help="The sequencing library insert size"
    )
    _parser.add_argument(
        "--sam-format",
        required=True,
        choices=["SAM", "BAM"],
        type=str,
        default="BAM",
        help="The sam file format, default: BAM"
    )

    _parser.add_argument(
        "--marker-matrix",
        required=True,
        help="IGC 9.9M marker matrix reference file."
    )

    _parser.add_argument(
        "--outdir",
        required=False,
        default="./",
        help="The directory for output of results."
    )

    _parser.add_argument(
        "--identity",
        required=False,
        type=float,
        default=0.95,
        help="The identity threshold used for filtering of result of alignment."
    )

    _parser.add_argument(
        "--conserved-list",
        required=False,
        default=None,
        help="List of conserved read that could be aligned to two or more references."
    )

    _parser.add_argument(
        "--suffix",
        required=False,
        default="",
        help="The suffix of output files. For rerunning of scripts with different parameters."
    )

    _parser.add_argument(
        "--output-type",
        required=False,
        choices=["evaluation", "standard", "non-zero", "all"],
        default="evaluation",
        help="The output file type of abundance calculation.\
        \"evaluation\": simplified format for abundance evaluation.\
        \"standard\": original abundance output similar to cOMG.\
        \"non-zero\": simplified standard output without zero value.\
        \"all\": output two files, one for evaluation and another for standard.\
        "
    )

    return _parser.parse_args()


if __name__ == '__main__':
    args = checkArgv()
    conservedReadInstance = ConservedRead(
        conservedReadListFile=args.conserved_list)
    markerMatrixRef = MarkerMatrix(args.marker_matrix)

    # time_rec1=time.time()
    sample = Sample(args.sample_name, args.sam_file, args.insert_size, args.sam_format,
                    markerMatrixRef, conservedReadInstance, identity=args.identity)
    # time_rec2=time.time()-time_rec1
    #sys.stderr.write("time for initialize sample:\t{0}s\n".format(time_rec2))

    # time_rec1=time.time()
    sample.ReadsCount()
    # time_rec2=time.time()-time_rec1
    #sys.stderr.write("time for readscount:\t{0}s\n".format(time_rec2))

    # time_rec1=time.time()
    sample.AbunCalculate()
    # time_rec2=time.time()-time_rec1
    #sys.stderr.write("time for abundance calculate:\t{0}s\n".format(time_rec2))

    # time_rec1=time.time()
    sample.AbundanceOutput(args.outdir, args.output_type, args.suffix)
    # time_rec2=time.time()-time_rec1
    #sys.stderr.write("time for abundance output:\t{0}s\n".format(time_rec2))
