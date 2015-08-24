
import sys
import getopt
import re
import copy


def usage():
    usagestr = """ python VCFnRefToFasta.py [options] -v data.vcf -r reference.fasta -o output_prefix
###############
-v | --vcf=             input vcf file
-r | --reference=       reference fasta file (presumably genomic)
-o | --out=             output_prefix (will be named .fasta, and if -g is set .transcript.fasta and .gff)
-q | --q_phred=         phred quality threshold (default = 20)
-g | --gff=             reference gff file (if supplied, script will return modified gff and associated transcriptome)
Note: -g is not yet implemented!
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


class VCFline:
    def __init__(self, vcfline_str):
        columns = vcfline_str.rstrip().split('\t')
        self.pos = int(columns[1])
        self.seqid = columns[0]
        self.ref = columns[3]
        self.alt = columns[4].split(',')
        self.info = columns[7]
        self.quality = self.get_quality()
        # note, the normal quality column is not used because it seems to be messed up in the vcf of interest
        # instead FQ from info

    def get_quality(self):
        pre_quality = re.search('FQ=([\-0-9]*)', self.info)
        # FQ for my gff of interest is the "Phred probability of all samples being the same"
        # nevermind that pred "probability" is not... the way I'd put it
        # appears to be the Phred score for variant call quality (which is reasonable)
        # anyways, you might have to write your own function for your appropriate quality score
        try:
            quality = int(pre_quality.group(1))
        except Exception:
            print pre_quality
            print self.info + " NO FQ FOUND"
        return quality


def split_seq(seq, l):
    """
    splits sequence (str) into list with subsequences of length (l)
    :param seq: list
    :param l: int
    :return: list of lists length l
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
    threshold = 20
    # get opt
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "v:r:g:o:q:h",
                                       ["vcf=", "reference=", "gff=", "out=", "q_phred=", "help"])
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
        elif o in ("-q", "--q_phred"):
            threshold = a
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if any([x is None for x in [vcfin, refin, fileout]]):
        print("input fasta required")
        usage()

    fastaout = open(fileout, 'w')
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
    for line in vcf.readlines(): #todo the right way 'with file as f...'
        if not line.startswith('#'):
            print line
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
                # update previous position to start after replacement
                previous_pos = end
                # keep track of total change in length for gff
                delta_l += len(v.alt[0]) - len(v.ref)
                print delta_l
            else:
                print '"' + fasta_dict[v.seqid][start:end] + '" where we expected: "' + v.ref + '"'
                raise Exception("mismatch to reference, REF not at POS in CHROM")

    for id in newfasta:
        fastaout.writelines('>' + id + '\n')
        for subseq in split_seq(newfasta[id], 60):
#            fastaout.writelines(''.join(subseq) + '\n') #todo break seq up by 60 and print all
            fastaout.writelines(subseq + '\n')

if __name__ == "__main__":
    main()
    test = list('abc')
    print test
    new = ['x', 'y']
    for item in reversed(new):
        test.insert(2, item)
#    test.insert(2,['x','y'])
    print test
