import os
import errno
import re

import edlib


'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        cigar_tuples.append((int(length), type_ ))

    r_index = 0
    q_index = 0
    q_aln = []
    r_aln = []
    for length_ , type_ in cigar_tuples:
        if type_ == "=" or type_ == "X":
            q_aln.append(query[q_index : q_index + length_])
            r_aln.append(ref[r_index : r_index + length_])

            r_index += length_
            q_index += length_
        
        elif  type_ == "I":
            # insertion w.r.t. reference
            r_aln.append('-' * length_)
            q_aln.append(query[q_index: q_index + length_])
            #  only query index change
            q_index += length_

        elif type_ == 'D':
            # deletion w.r.t. reference
            r_aln.append(ref[r_index: r_index + length_])
            q_aln.append('-' * length_)
            #  only ref index change
            r_index += length_
        
        else:
            print("error")
            print(cigar)
            sys.exit()

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln])

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def read_barcodes(primer_file):
    barcodes = {}
    for i, line in enumerate(open(primer_file, "r")):
        if i == 0:
            continue
        # acc, FwIndex, FWPrimer, RvIndex, RvPrimer = line.split()
        # barcodes[acc+"_FwIndex_fwd"] = FwIndex.upper()
        # barcodes[acc+ "_RvIndex_fwd"] = RvIndex.upper()

        # FwIndex_rev = reverse_complement(FwIndex.upper())
        # RvIndex_rev = reverse_complement(RvIndex.upper())
        # barcodes[acc+"_FwIndex_rev"] = FwIndex_rev
        # barcodes[acc+ "_RvIndex_rev"] = RvIndex_rev

        acc, FwIndex, FWPrimer, RvIndex, RvPrimer = line.split()
        barcodes[acc+"_FwIndex_fwd"] = FWPrimer.upper()[:25]
        barcodes[acc+ "_RvIndex_fwd"] = RvPrimer.upper()[:25]

        FwIndex_rev = reverse_complement(FWPrimer.upper()[:25])
        RvIndex_rev = reverse_complement(RvPrimer.upper()[:25])
        barcodes[acc+"_FwIndex_rev"] = FwIndex_rev
        barcodes[acc+ "_RvIndex_rev"] = RvIndex_rev

    return barcodes

def get_universal_tails():
    barcodes = {'1_fwd' : 'TTTCTGTTGGTGCTGATATTGC',
                 '2_fwd' : 'ACTTGCCTGTCGCTCTATCTTC'}
    barcodes['1_rev'] = reverse_complement(barcodes['1_fwd'])
    barcodes['2_rev'] = reverse_complement(barcodes['2_fwd'])
    return barcodes


def find_barcode_locations(center, barcodes, primer_max_ed):
    "Find barcodes in a center using edlib"
    all_locations = []
    for primer_acc, primer_seq in barcodes.items():
        # print(primer_acc, primer_seq,center)
        result = edlib.align(primer_seq, center,
                             mode="HW", task="locations", k=primer_max_ed)
        ed = result["editDistance"]
        locations = result["locations"]
        # print(locations, ed)
        if locations:
            all_locations.append((primer_acc, locations[0][0], locations[0][1], ed))
    #     if locations:
    #         # all_locations[primer_acc] = []
    #         for refstart, refend in locations:
    #             refend += 1
    #             # ('Hit', 'Ref RefStart RefEnd Query QueryStart QueryEnd Score')
    #             hit = Hit(read.Name, refstart, refend, primer_acc, 0,
    #                       len(primer_seq),  ed / len(primer_seq))
    #             all_locations.append(hit)
    # refined_locations = refine_locations(read, all_primers, all_locations)
    return all_locations







