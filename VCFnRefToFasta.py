
import sys
import getopt
import re
import copy
import intervaltree


def usage():
    usagestr = """ python VCFnRefToFasta.py [options] -v data.vcf -r reference.fasta -o output_prefix
###############
-v | --vcf=             input vcf file
-r | --reference=       reference fasta file (presumably genomic)
-o | --out=             output_prefix (will be named .fa, and if -g is set and .gff)
-q | --q_phred=         phred quality threshold (default = 20)
-g | --gff=             reference gff file (if supplied, script will return modified gff and associated transcriptomei)
-h | --help             prints this message

Note this is probably quite specific for the VCF produced in the 150 tomato genomes project. But that's what I needed.
"""
    print(usagestr)
    sys.exit(1)


def fasta_to_dict(fasta_str):
    fasta_dict = {}
    last_id = None
    seq = ''
    for line in fasta_str.split('\n'):
        if line.startswith('>'):
            if last_id is not None:
                fasta_dict[last_id] = seq
            last_id = line.replace('>', '')
            seq = ''
        elif line != '':
            seq += line
#    fasta_dict[last_id] = list(seq)
    fasta_dict[last_id] = seq
    return fasta_dict


def lists_match(list1, list2):
    if len(list1) != len(list2):
        out = False
    else:
        out = all([list1[i] == list2[i] for i in range(len(list1))])
    return out


class Offsets:
    def __init__(self, status='matching', offset=0):
        if status in ['matching', 'not-matching']:
            self.status = status
        else:
            raise ValueError('status accepts "matching" or "not-matching" only')
        self.offset = offset


class Gffentry:
    def __init__(self, line):
        splitline = line.rstrip().split('\t')
        self.seqname = splitline[0]
        self.source = splitline[1] + ',vcfhd'
        self.feature = splitline[2]
        self.start = int(splitline[3])
        self.end = int(splitline[4])
        self.score = splitline[5]
        self.strand = splitline[6]
        self.frame = splitline[7]
        self.attribute = splitline[8]

    def __str__(self):
        splitline = [self.seqname, self.source, self.feature, self.start, self.end, self.score, self.strand,
                     self.frame, self.attribute]
        splitline = [str(x) for x in splitline]
        out = '\t'.join(splitline)
        return out


class VCFline:
    def __init__(self, vcfline_str):
        columns = vcfline_str.rstrip().split('\t')
        self.pos = int(columns[1])
        self.seqid = columns[0]
        self.ref = columns[3]
        self.alt = columns[4].split(',')
        self.info = columns[7]
        self.quality = float(columns[5])


def split_seq(seq, l):
    """
    splits sequence (str) into list with subsequences of length (l)
    :param seq: list
    :param l: int
    :return: list of strs length l
    """
    out = []
    for i in range(0, len(seq), l):
        upper = min(i + l, len(seq))
        subseq = seq[i:upper]
        out += [subseq]
    return out


def main():
    vcfin = None
    refin = None
    gffin = None
    fileout = None
    threshold = 10
    # get opt
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "v:r:g:o:q:h",
                                       ["vcf=", "reference=", "gff=", "out=", "quality=", "help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()

    for o, a in opts:
        if o in ("-v", "--vcf"):
            vcfin = a
        elif o in ("-r", "--reference"):
            refin = a
        elif o in ("-g", "--gff"):
            gffin = a
        elif o in ("-o", "--out"):
            fileout = a
        elif o in ("-q", "--quality"):
            threshold = a
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if any([x is None for x in [vcfin, refin, fileout]]):
        print("input vcf, reference fasta, and output prefix required")
        usage()

    fastaout = open(fileout + '.fa', 'w')
    vcf = open(vcfin)
    fasta = open(refin)

    fasta_dict = fasta_to_dict(fasta.read().rstrip())
    print fasta_dict.keys()
    newfasta = {}
    for key in fasta_dict:
        newfasta[key] = ''#list(fasta_dict[key])
#    newfasta = copy.deepcopy(fasta_dict)
    previous_pos = 0
    delta_l = 0
    offset_tree = intervaltree.IntervalTree()
    for line in vcf.readlines(): #todo the right way 'with file as f...'
        if not line.startswith('#'):
#            print line
            v = VCFline(line)
            reflist = list(v.ref)
            # check position (-1 because vcf counts from 1 not 0)
            start = v.pos - 1
            end = start + len(v.ref)
            if fasta_dict[v.seqid][start:end] == v.ref:
                # add from last change to just before replacement
                newfasta[v.seqid] += fasta_dict[v.seqid][previous_pos:v.pos]
                # add alternative
                newfasta[v.seqid] += v.alt[0]
                # todo, something about heterozygosity and cases where you don't simply want the first listed
                # tracking for gff
                between = Offsets('matching', delta_l)
                overlap = Offsets('not-matching', delta_l)
                offset_tree[previous_pos:start] = between
                offset_tree[start:end] = overlap
                # update previous position to start after replacement
                previous_pos = end
                # keep track of total change in length for gff
                delta_l += len(v.alt[0]) - len(v.ref)
#                print delta_l
            else:
                print '"' + fasta_dict[v.seqid][start:end] + '" where we expected: "' + v.ref + '"'
                raise Exception("mismatch to reference, REF not at POS in CHROM")

    if gffin is not None:
        gfflines = open(gffin).readlines()
        gffout = open(fileout + '.gff', 'w')
        for line in gfflines:
            line = line.rstrip()
            if not line.startswith('##'):
                gffentry = Gffentry(line)
                offsetstart = sorted(offset_tree[gffentry.start])[0]
                gffentry.start += offsetstart.data.offset
                offsetend = sorted(offset_tree[gffentry.end])[0]
                gffentry.end += offsetend.data.offset
                gffout.writelines(str(gffentry) + '\n')
            elif line.find('annot-version') > -1:
                gffout.writelines(line + '-modified_with_vcfhd' + '\n')  # todo, better annotation version update
            else:
                gffout.writelines(line + '\n')

    for seqid in newfasta:
        fastaout.writelines('>' + seqid + '\n')
        for subseq in split_seq(newfasta[seqid], 60):
            fastaout.writelines(subseq + '\n')

if __name__ == "__main__":
    main()
