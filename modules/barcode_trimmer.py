
import edlib

from modules import help_functions

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def read_barcodes(primer_file):
    barcodes = { acc + '_fw' : seq.strip() for acc, (seq, _) in help_functions.readfq(open(primer_file, 'r'))}

    for acc, seq in list(barcodes.items()):
        print(acc, seq,acc[:-3])
        barcodes[acc[:-3] + '_rc'] = reverse_complement(seq.upper())

    print(barcodes)
    return barcodes

def get_universal_tails():
    barcodes = {'1_F_fw' : 'TTTCTGTTGGTGCTGATATTGC',
                 '2_R_rc' : 'ACTTGCCTGTCGCTCTATCTTC'}
    barcodes['1_F_rc'] = reverse_complement(barcodes['1_F_fw'])
    barcodes['2_R_fw'] = reverse_complement(barcodes['2_R_rc'])
    print(barcodes)
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
        print(locations, ed)
        if locations:
            all_locations.append((primer_acc, locations[0][0], locations[0][1], ed))
    return all_locations


def remove_barcodes(centers, barcodes, args):
    """
        Modifies consensus sequences by copping of at barcode sites.
        This implies changing the datastructure centers with the modified consensus sequeces
    """
    for i, (nr_reads_in_cluster, c_id, center, reads_path_name) in enumerate(centers):
        barcode_locations = find_barcode_locations(center, barcodes, args.primer_max_ed) 
        if barcode_locations:
            print("FOUND BARCODE", barcode_locations)
            cut_start = 0
            cut_end = len(center)
            print(center)
            for bc, start, stop, ed in barcode_locations:
                # print(ed,bc, bc[-4], bc[-2:])
                if bc[-4] == 'F' and bc[-2:] == 'fw':
                    cut_start = stop
                elif bc[-4] == 'R' and bc[-2:] == 'fw': 
                    cut_end = start
                elif bc[-4] == 'R' and bc[-2:] == 'rc':
                    cut_start = stop
                elif bc[-4] == 'F' and bc[-2:] == 'rc': 
                    cut_end = start
                else:
                    print()
                    print("Primer file not in correct format!")
                    print()
            # print(center)
            center = center[cut_start: cut_end]
            print(center, "NEW")
            print("cut start", cut_start, "cut end", cut_end)
            centers[i][2] = center
