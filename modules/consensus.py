from __future__ import print_function
import subprocess
import sys, os
from sys import stdout
import re 
import shutil
import parasail
import glob

from modules import help_functions


def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    cig_pos = 0
    for length in result[:-1]:
        cig_pos += len(length)
        type_ = cigar[cig_pos]
        cig_pos += 1
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

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln]), cigar_tuples


def parasail_alignment(s1, s2, match_score = 2, mismatch_penalty = -2, opening_penalty = 3, gap_ext = 1):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(s1, s2, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        print("SATURATED!",len(s1), len(s2))
        result = parasail.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, user_matrix)
        print("computed 32 bit instead")

    # difference in how to obtain string from parasail between python v2 and v3... 
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    s1_alignment, s2_alignment, cigar_tuples = cigar_to_seq(cigar_string, s1, s2)
    # print(result.score, len(s1), len(s2))
    # print(s1_alignment)
    # print(s2_alignment)
    # print(cigar_string)
    # sys.exit()

    return s1_alignment, s2_alignment, cigar_string, cigar_tuples, result.score

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def run_spoa(reads, spoa_out_file, spoa_path):
    with open(spoa_out_file, "w") as output_file:
        # print('Running spoa...', end=' ')
        stdout.flush()
        with open("/dev/null", "w") as null:
            subprocess.check_call([ spoa_path, reads, "-l", "0", "-r", "0", "-g", "-2"], stdout=output_file, stderr=null)
        # print('Done.')
        stdout.flush()
    with open(spoa_out_file, "r") as sof:
        l = sof.readlines()
    consensus = l[1].strip()
    return consensus

def run_medaka(reads_to_center, center_file, outfolder, cores, medaka_model):
    medaka_stdout = os.path.join(outfolder, "stdout.txt")
    with open(medaka_stdout, "w") as output_file:
        # print('Running medaka...', end=' ')
        stdout.flush()
        with open(os.path.join(outfolder, "stderr.txt"), "w") as medaka_stderr:
            if medaka_model:
                subprocess.check_call(['medaka_consensus', '-i', reads_to_center, "-d", center_file, "-o", outfolder, "-t", cores, "-m", medaka_model], stdout=output_file, stderr=medaka_stderr)
            else:
                subprocess.check_call(['medaka_consensus', '-i', reads_to_center, "-d", center_file, "-o", outfolder, "-t", cores], stdout=output_file, stderr=medaka_stderr)

        # print('Done.')
        stdout.flush()

def run_racon(reads_to_center, center_file, outfolder, cores, racon_iter):
    racon_stdout = os.path.join(outfolder, "stdout.txt")
    with open(racon_stdout, "w") as output_file:
        # print('Running medaka...', end=' ')
        stdout.flush()
        for i in range(racon_iter):
            with open(
                os.path.join(outfolder, "read_alignments_it_{0}.paf".format(i)), 'w'
            ) as read_alignments, open(
                os.path.join(outfolder, "mm2_stderr_it_{0}.txt".format(i)), "w"
            ) as mm2_stderr, open(
                os.path.join(outfolder, "racon_stderr_it_{0}.txt".format(i)), "w"
            ) as racon_stderr, open(
                os.path.join(outfolder, "racon_polished_it_{0}.fasta".format(i)
            ), 'w') as racon_polished:
                subprocess.check_call(['minimap2', '-x', 'map-ont', center_file, reads_to_center], stdout=read_alignments, stderr=mm2_stderr)
                subprocess.check_call(['racon', reads_to_center, read_alignments.name, center_file], stdout=racon_polished, stderr=racon_stderr)
            center_file = racon_polished.name

        shutil.copyfile(center_file, os.path.join(outfolder, "consensus.fasta"))
        # print('Done.')
        stdout.flush()
    # consensus.run_medaka( " medaka_consensus -i ~/tmp/stefan_isonclust/mixed_b1_b3_b4_b5.fastq -d ~/tmp/stefan_isonclust/mixed_b1_b3_b4_b5_isonclust/consensus_references.fasta -o ~/tmp/stefan_isonclust/mixed_b1_b3_b4_b5_medaka/ -t 1 -m r941_min_high_g303")


def highest_aln_identity(seq, seq2):
    # RC
    seq2_rc = reverse_complement(seq2)
    seq_aln_rc, seq2_aln_rc, cigar_string_rc, cigar_tuples_rc, alignment_score_rc = parasail_alignment(seq, seq2_rc)
    nr_mismatching_pos = len([1 for n1, n2 in zip(seq_aln_rc, seq2_aln_rc) if n1 != n2])
    total_pos_rc = len(seq_aln_rc)
    aln_identity_rc =  (total_pos_rc - nr_mismatching_pos) / float(total_pos_rc)
    print('Rec comp orientation identity %: ', aln_identity_rc)

    # FW
    seq_aln, seq2_aln, cigar_string, cigar_tuples, alignment_score = parasail_alignment(seq, seq2)
    nr_mismatching_pos = len([1 for n1, n2 in zip(seq_aln, seq2_aln) if n1 != n2])
    total_pos = len(seq_aln)
    aln_identity_fw = (total_pos - nr_mismatching_pos) / float(total_pos)  
    print('Forward orientation identity %: ', aln_identity_fw)
    aln_identity = max([aln_identity_fw, aln_identity_rc])
    return aln_identity


def detect_reverse_complements(centers, rc_identity_threshold):
    filtered_centers = []
    already_removed = set()
    for i, (nr_reads_in_cl, c_id, seq, reads_path) in enumerate(centers):
        if type(reads_path) != list:
            all_reads = [reads_path]
        else:
            all_reads = reads_path

        merged_cluster_id = c_id
        merged_nr_reads = nr_reads_in_cl
        if c_id in already_removed:
            print("has already been merged, skipping")
            continue

        elif i == len(centers) - 1: # last sequence and it is not in already_removed
            filtered_centers.append( [merged_nr_reads, c_id, seq, all_reads ] )

        else:
            for j, (nr_reads_in_cl2, c_id2, seq2, reads_path) in enumerate(centers[i+1:]):
                aln_identity = highest_aln_identity(seq, seq2)
                if aln_identity >= rc_identity_threshold:
                    print("Detected two consensus sequences with alignment identidy above threshold (from either reverse complement or split clusters). Keeping center with the most read support and merging reads.")
                    merged_nr_reads += nr_reads_in_cl2
                    already_removed.add(c_id2)

                    if type(reads_path) != list:
                        all_reads.append(reads_path)
                    else:
                        for rp in reads_path:
                            all_reads.append(rp)

            filtered_centers.append( [merged_nr_reads, c_id, seq, all_reads] )

    print(len(filtered_centers), "consensus formed.")
    return filtered_centers


def polish_sequences(centers, args):
    print("Saving spoa references to files:", os.path.join(args.outfolder, "consensus_reference_X.fasta"))
    # printing output from spoa and grouping reads
    # to_polishing = []
    if args.medaka:
        polishing_pattern = os.path.join(args.outfolder, "medaka_cl_id_*")
    elif args.racon:
        polishing_pattern = os.path.join(args.outfolder, "racon_cl_id_*")

    for folder in glob.glob(polishing_pattern):
        shutil.rmtree(folder)

    spoa_pattern = os.path.join(args.outfolder, "consensus_reference_*")
    for file in glob.glob(spoa_pattern):
        os.remove(file)

    for i, (nr_reads_in_cluster, c_id, center, all_reads) in enumerate(centers):
        # print('lol',c_id,center)
        spoa_center_file = os.path.join(args.outfolder, "consensus_reference_{0}.fasta".format(c_id))
        with open(spoa_center_file, "w") as f:
            f.write(">{0}\n{1}\n".format("consensus_cl_id_{0}_total_supporting_reads_{1}".format(c_id, nr_reads_in_cluster), center))
        
        nr_reads_used = 0
        all_reads_file = os.path.join(args.outfolder, "reads_to_consensus_{0}.fastq".format(c_id))
        with open(all_reads_file, "w") as f:
            for fasta_file in all_reads: 
                reads = { acc : (seq, qual) for acc, (seq, qual) in help_functions.readfq(open(fasta_file, 'r'))}
                for acc, (seq, qual) in reads.items():
                    acc_tmp = acc.split()[0]
                    f.write("@{0}\n{1}\n{2}\n{3}\n".format(acc_tmp, seq, "+", qual))
                    nr_reads_used += 1
        # to_polishing.append( (nr_reads_in_cluster, c_id, spoa_center_file, all_reads_file) )

        if args.medaka:
            print("running medaka on spoa reference {0} using {1} reads for polishing.".format(c_id, nr_reads_used))
            # for (nr_reads_in_cluster, c_id, spoa_center_file, all_reads_file) in to_polishing:
            polishing_outfolder = os.path.join(args.outfolder, "medaka_cl_id_{0}".format(c_id))
            help_functions.mkdir_p(polishing_outfolder)
            run_medaka(all_reads_file, spoa_center_file, polishing_outfolder, "1", args.medaka_model)
            print("Saving medaka reference to file:", os.path.join(args.outfolder, "medaka_cl_id_{0}/consensus.fasta".format(c_id)))   
            with open(os.path.join(polishing_outfolder, "consensus.fasta"), 'r') as cf:
                l = cf.readlines()
            center_polished = l[1].strip()
            centers[i][2] = center_polished
        elif args.racon:
            print("running racon on spoa reference {0} using {1} reads for polishing.".format(c_id, nr_reads_used))
            # for (nr_reads_in_cluster, c_id, spoa_center_file, all_reads_file) in to_polishing:
            polishing_outfolder = os.path.join(args.outfolder, "racon_cl_id_{0}".format(c_id))
            help_functions.mkdir_p(polishing_outfolder)
            run_racon(all_reads_file, spoa_center_file, polishing_outfolder, "1", args.racon_iter)
            print("Saving racon reference to file:", os.path.join(args.outfolder, "racon_cl_id_{0}/consensus.fasta".format(c_id)))   
            with open(os.path.join(polishing_outfolder, "consensus.fasta"), 'r') as cf:
                l = cf.readlines()
            center_polished = l[1].strip()
            centers[i][2] = center_polished

    return centers


def form_draft_consensus(clusters, representatives, sorted_reads_fastq_file, work_dir, abundance_cutoff, args):
    centers = []
    singletons = 0
    discarded_clusters = []
    reads = { acc : (seq, qual) for acc, (seq, qual) in help_functions.readfq(open(sorted_reads_fastq_file, 'r'))}
    for c_id, all_read_acc in sorted(clusters.items(), key = lambda x: (len(x[1]),representatives[x[0]][5]), reverse=True):
        nr_reads_in_cluster = len(all_read_acc)
        # print("nr_reads_in_cluster", nr_reads_in_cluster)
        if nr_reads_in_cluster >= abundance_cutoff:
            reads_path_name = os.path.join(work_dir, "reads_c_id_{0}.fq".format(c_id))
            with open(reads_path_name, "w") as reads_file:
                for i, acc in enumerate(all_read_acc):
                    if (args.max_seqs_for_consensus) >=0 and (i >= args.max_seqs_for_consensus):
                        break
                    seq, qual = reads[acc]
                    reads_file.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, seq, "+", qual))
                    # reads.write(">{0}\n{1}\n".format(str(q_id)+str(pos1)+str(pos2), seq))
            # spoa_ref = create_augmented_reference.run_spoa(reads_path_name, os.path.join(work_dir,"spoa_tmp.fa"), "spoa")
            tmp_param = args.max_seqs_for_consensus if args.max_seqs_for_consensus > 0 else 2**32
            print("creating center of {0} sequences.".format(min(nr_reads_in_cluster, tmp_param)))
            center = run_spoa(reads_path_name, os.path.join(work_dir,"spoa_tmp.fa"), "spoa")
            centers.append( [nr_reads_in_cluster, c_id, center, reads_path_name])
        elif nr_reads_in_cluster == 1:
            singletons += 1
        elif nr_reads_in_cluster > 1:
            discarded_clusters.append(nr_reads_in_cluster)
    print(f"{singletons} singletons were discarded")
    print(
        f"{len(discarded_clusters)} clusters were discarded due to not passing the abundance_cutoff: "
        f"a total of {sum(discarded_clusters)} reads were discarded. "
        f"Highest abundance among them: {max(discarded_clusters or [0])} reads."
    )
    return centers
