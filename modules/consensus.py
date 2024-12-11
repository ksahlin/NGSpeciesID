from __future__ import print_function
import subprocess
import sys, os
from sys import stdout
import re 
import shutil
import parasail
import glob
import logging

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
            logging.error("Error processing cigar")
            logging.error(cigar)
            sys.exit()

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln]), cigar_tuples


def parasail_alignment(s1, s2, match_score = 2, mismatch_penalty = -2, opening_penalty = 3, gap_ext = 1):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(s1, s2, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        logging.warning(f"SATURATED!{len(s1)} {len(s2)}")
        result = parasail.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, user_matrix)
        logging.warning("computed 32 bit instead")

    # difference in how to obtain string from parasail between python v2 and v3... 
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    s1_alignment, s2_alignment, cigar_tuples = cigar_to_seq(cigar_string, s1, s2)

    return s1_alignment, s2_alignment, cigar_string, cigar_tuples, result.score

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def run_spoa(reads, spoa_out_file, spoa_path):
    with open(spoa_out_file, "w") as output_file:
        stdout.flush()
        with open("/dev/null", "w") as null:
            subprocess.check_call([ spoa_path, reads, "-l", "0", "-r", "0", "-g", "-2"], stdout=output_file, stderr=null)
        stdout.flush()
    with open(spoa_out_file, "r") as sof:
        l = sof.readlines()
    consensus = l[1].strip()
    return consensus

def run_medaka(reads_to_center, center_file, outfolder, cores, medaka_model, outfastq=False):
    medaka_stdout = os.path.join(outfolder, "stdout.txt")
    with open(medaka_stdout, "w") as output_file:
        stdout.flush()
        with open(os.path.join(outfolder, "stderr.txt"), "w") as medaka_stderr:
            cmd_args = ['medaka_consensus', '-i', reads_to_center, "-d", center_file, "-o", outfolder, "-t", cores]
            if medaka_model:
                cmd_args += ["-m", medaka_model]
            if outfastq:
                cmd_args += ["-q"]
            subprocess.check_call(cmd_args, stdout=output_file, stderr=medaka_stderr)
        stdout.flush()

def run_racon(reads_to_center, center_file, outfolder, cores, racon_iter):
    racon_stdout = os.path.join(outfolder, "stdout.txt")
    with open(racon_stdout, "w") as output_file:
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
        stdout.flush()


def highest_aln_identity(seq, seq2):
    # RC
    seq2_rc = reverse_complement(seq2)
    seq_aln_rc, seq2_aln_rc, cigar_string_rc, cigar_tuples_rc, alignment_score_rc = parasail_alignment(seq, seq2_rc)
    nr_mismatching_pos = len([1 for n1, n2 in zip(seq_aln_rc, seq2_aln_rc) if n1 != n2])
    total_pos_rc = len(seq_aln_rc)
    aln_identity_rc =  (total_pos_rc - nr_mismatching_pos) / float(total_pos_rc)
    logging.debug(f"Rec comp orientation identity %: {aln_identity_rc}")

    # FW
    seq_aln, seq2_aln, cigar_string, cigar_tuples, alignment_score = parasail_alignment(seq, seq2)
    nr_mismatching_pos = len([1 for n1, n2 in zip(seq_aln, seq2_aln) if n1 != n2])
    total_pos = len(seq_aln)
    aln_identity_fw = (total_pos - nr_mismatching_pos) / float(total_pos)  
    logging.debug(f"Forward orientation identity %: {aln_identity_fw}")
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
            logging.debug("has already been merged, skipping")
            continue

        elif i == len(centers) - 1: # last sequence and it is not in already_removed
            filtered_centers.append( [merged_nr_reads, c_id, seq, all_reads ] )

        else:
            for j, (nr_reads_in_cl2, c_id2, seq2, reads_path) in enumerate(centers[i+1:]):
                aln_identity = highest_aln_identity(seq, seq2)
                if aln_identity >= rc_identity_threshold:
                    logging.debug("Detected two consensus sequences with alignment identidy above threshold (from either reverse complement or split clusters). Keeping center with the most read support and merging reads.")
                    merged_nr_reads += nr_reads_in_cl2
                    already_removed.add(c_id2)

                    if type(reads_path) != list:
                        all_reads.append(reads_path)
                    else:
                        for rp in reads_path:
                            all_reads.append(rp)

            filtered_centers.append( [merged_nr_reads, c_id, seq, all_reads] )

    logging.debug(f"{len(filtered_centers)} consensus formed.")
    return filtered_centers


def polish_sequences(centers, args):
    spoa_ref_location = os.path.join(args.outfolder, "consensus_reference_X.fasta")
    logging.debug(f"Saving spoa references to files: {spoa_ref_location}")
    # printing output from spoa and grouping reads
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

        if args.medaka:
            logging.debug("running medaka on spoa reference {0} using {1} reads for polishing.".format(c_id, nr_reads_used))
            polishing_outfolder = os.path.join(args.outfolder, "medaka_cl_id_{0}".format(c_id))
            outfiles = [  # consider all output formats for compatibility with all Medaka versions
                os.path.join(polishing_outfolder, "consensus.fasta"),
                os.path.join(polishing_outfolder, "consensus.fastq")
            ]
            help_functions.mkdir_p(polishing_outfolder)
            run_medaka(all_reads_file, spoa_center_file, polishing_outfolder, "1", args.medaka_model, outfastq=args.medaka_fastq)
            medaka_ref_location = os.path.join(polishing_outfolder, "consensus.fasta/q")
            logging.debug(f"Saving medaka reference to file: {medaka_ref_location}")
            for f in outfiles:
                if os.path.isfile(f):
                    with open(f, 'r') as cf:
                        centers[i][2] = cf.readlines()[1].strip()  # the second line is the nucleotide sequence
                        break
            assert centers[i][2], "Medaka consensus sequence not found"
        elif args.racon:
            logging.debug("running racon on spoa reference {0} using {1} reads for polishing.".format(c_id, nr_reads_used))
            polishing_outfolder = os.path.join(args.outfolder, "racon_cl_id_{0}".format(c_id))
            help_functions.mkdir_p(polishing_outfolder)
            run_racon(all_reads_file, spoa_center_file, polishing_outfolder, "1", args.racon_iter)
            racon_ref_location = os.path.join(args.outfolder, "racon_cl_id_{0}/consensus.fasta".format(c_id))
            logging.debug(f"Saving racon reference to file: {racon_ref_location}")
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
        if nr_reads_in_cluster >= abundance_cutoff:
            reads_path_name = os.path.join(work_dir, "reads_c_id_{0}.fq".format(c_id))
            with open(reads_path_name, "w") as reads_file:
                for i, acc in enumerate(all_read_acc):
                    if (args.max_seqs_for_consensus) >=0 and (i >= args.max_seqs_for_consensus):
                        break
                    seq, qual = reads[acc]
                    reads_file.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, seq, "+", qual))
            tmp_param = args.max_seqs_for_consensus if args.max_seqs_for_consensus > 0 else 2**32
            logging.debug("creating center of {0} sequences.".format(min(nr_reads_in_cluster, tmp_param)))
            center = run_spoa(reads_path_name, os.path.join(work_dir,"spoa_tmp.fa"), "spoa")
            centers.append( [nr_reads_in_cluster, c_id, center, reads_path_name])
        elif nr_reads_in_cluster == 1:
            singletons += 1
        elif nr_reads_in_cluster > 1:
            discarded_clusters.append(nr_reads_in_cluster)
    logging.debug(f"{singletons} singletons were discarded")
    logging.debug(
        f"{len(discarded_clusters)} clusters were discarded due to not passing the abundance_cutoff: "
        f"a total of {sum(discarded_clusters)} reads were discarded. "
        f"Highest abundance among them: {max(discarded_clusters or [0])} reads."
    )
    return centers
