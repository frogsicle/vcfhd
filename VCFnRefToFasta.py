
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
        elif line != '':
            seq += line
    return fasta_dict


class VCFline:
    def __init__(self, vcfline_str):
        columns = vcfline_str.rstrip().split('\t')
        self.pos = columns[1]
        self.seqid = columns[0]
        self.ref = columns[3]
        self.alt = columns[4]
        self.info = columns[0]
        self.quality = self.get_quality(self)
        # note, the normal quality column is not used because it seems to be messed up in the vcf of interest
        # instead FQ from info

    def get_quailty(self):
        pre_quality = re.search('FQ=([\-0-9]*);]')
        # FQ for my gff of interest is the "Phred probability of all samples being the same"
        # nevermind that pred "probability" is not... the way I'd put it
        # appears to be the Phred score for variant call quality (which is reasonable)
        # anyways, you might have to write your own function for your appropriate quality score
        quality = int(pre_quality)
        return quality


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
    newfasta = copy.deepcopy(fasta_dict)
    for line in reversed(vcf.readlines()):  # todo reversed efficient line by line opening of file
        v = VCFline(line)
        # check position
        if fasta_dict[v.seqid][v.pos:(v.pos + len(v.ref))] == v.ref:
            newfasta[v.seqid][v.pos:] = newfasta[v.seqid][v.pos].replace(v.ref, v.alt, 1)
        else:
            raise Exception("mismatch to reference, REF not at POS in CHROM, line: \n" + line)
    for id in newfasta:
        fastaout.writelines('>' + id)
        fastaout.writelines(newfasta[id][:60]) #todo break seq up by 60 and print all

if __name__ == "__main__":
    main()
