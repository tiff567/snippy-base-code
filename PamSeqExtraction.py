# Pam gRNA sequence extraction based on Position from given Fasta data

from Bio.Seq import Seq

# spcas9 type for A->G mutation. Pam sequence: NGG
def get_pam_sequence_ag_spcas9(fasta_data, chr_num, chr_pos, shift_pos = 14, seq_len = 20, is_optimal = True):

    pam_num = 0         # count the # of pam sequences
    optimal_num = 0     # count the # of optimal pam sequences
    grna_str = ""       # gRNA string for stored Pam Sequences with space separation
    grna_pos_str = ""     # Pam Sequence Postion for stored Pam Sequence positions with space separation

    pam_checking_times = 4 # the total number of pam searching window shifts
    pam_checking_pos = shift_pos + pam_checking_times  # the max shift position

    while shift_pos < pam_checking_pos:
        isPam = False
        pam_str = fasta_data['chr' + str(chr_num)][chr_pos + shift_pos:chr_pos + shift_pos + 2]   # get 2 characters downstream
        if str(pam_str) == 'GG':  # check whether it is PAM
            isPam = True

        if isPam:
            # Get the PAM sequence from the PAM position upstream 20 characters
            grna = fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - seq_len:chr_pos + shift_pos]

            # check whether it is optimal
            if is_optimal:
                  if str(fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 4:chr_pos + shift_pos - 3]) == 'A' or str(
                        fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 4:chr_pos + shift_pos - 3]) == 'T':

                    grna_str += (str(grna) + " ")  # write done the gRNA string
                    grna_pos_str += (str(chr_pos + shift_pos - seq_len) + " ")  # write down PAM sequence position
                    optimal_num += 1
            else:
                grna_str += (str(grna) + " ")    # write done the gRNA string
                grna_pos_str += (str(chr_pos + shift_pos -seq_len) + " ")  # write down PAM sequence position

            pam_num += 1  # count the pam num
        shift_pos += 1

    return grna_str, grna_pos_str, pam_num, optimal_num

# spcas9 type for C->G mutation. Pam sequence: NCC
def get_pam_sequence_tc_spcas9(fasta_data, chr_num, chr_pos, shift_pos = 14, seq_len = 20, is_optimal = True):

    pam_num = 0  # count the # of pam sequences
    optimal_num = 0  # count the # of optimal pam sequences
    grna_str = ""  # gRNA string for stored Pam Sequences with space separation
    grna_pos_str = ""     # Pam Sequence Postion for stored Pam Sequence positions with space separation

    pam_checking_times = 4  # the total number of pam searching window shifts
    pam_checking_pos = shift_pos + pam_checking_times  # the max shift position

    while shift_pos < pam_checking_pos:
        isPam = False
        pam_str = fasta_data['chr' + str(chr_num)][chr_pos - shift_pos:chr_pos - shift_pos + 2]   # get 2 characters upstream

        if str(pam_str) == 'CC':
            isPam = True

        if isPam:

            # Get the PAM sequence from the PAM position upstream 20 characters
            grna = fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - seq_len:chr_pos - shift_pos]

            # check whether it is optimal
            if is_optimal:
                if str(fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - 4:chr_pos - shift_pos - 3]) == 'A' or \
                        str(fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - 4:chr_pos - shift_pos - 3]) == 'T':
                    grna_seq = Seq(str(grna))  # convert it to BioSequence
                    grna_tc = grna_seq.reverse_complement()
                    grna_str += (str(grna_tc) + " ")  # write done the gRNA sequence
                    grna_pos_str += (str(chr_pos - shift_pos -seq_len) + " ")  # write down PAM sequence position
                    optimal_num += 1
            else:
                grna_seq = Seq(str(grna))  # convert it to BioSequence
                grna_tc = grna_seq.reverse_complement()
                grna_str += (str(grna_tc) + " ")
                grna_pos_str += (str(chr_pos - shift_pos - seq_len) + " ")  # write down PAM sequence position

            pam_num += 1  # count the pam num
        shift_pos += 1

    return grna_str, grna_pos_str, pam_num, optimal_num

# sacas9 type for A-G mutation. Pam sequence: NGRRT or NGRRN  R = A or G
def get_pam_sequence_ag_sacas9(fasta_data, chr_num, chr_pos, shift_pos = 14, seq_len = 20, is_optimal = True):

    pam_num = 0  # count the # of pam sequences
    optimal_num = 0  # count the # of optimal pam sequences
    grna_str = ""  # gRNA string for stored Pam Sequences with space separation
    grna_pos_str = ""     # Pam Sequence Postion for stored Pam Sequence positions with space separation

    pam_checking_times = 4  # the total number of pam searching window shifts
    pam_checking_pos = shift_pos + pam_checking_times  # the max shift position

    while shift_pos < pam_checking_pos:
        isPam = False
        pam_str = fasta_data['chr' + str(chr_num)][chr_pos + shift_pos:chr_pos + shift_pos + 4]   # get 4 characters downsteam

        # check whether it is PAM NGRRT (NGAAT, NGAGT, NGGAT, NGGGT, NGAATN, NGAGN, NGGAN, NGGGN
        if (str(pam_str) == 'GAAT') or (str(pam_str) == 'GAGT') or (str(pam_str) == 'GGAT') or (str(pam_str) == 'GGGT') \
                or (str(pam_str[0:3]) == 'GAA') or (str(pam_str[0:3]) == 'GAG') or (str(pam_str[0:3]) == 'GGA') or (str(pam_str[0:3]) == 'GGG'):
            isPam = True

        if isPam:
            # Get the PAM sequence from the PAM position upstream 20 characters
            grna = fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - seq_len:chr_pos + shift_pos]

            # check whether it is optimal
            if is_optimal:
                if str(fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 4:chr_pos + shift_pos - 3]) == 'A' or str(
                        fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 4:chr_pos + shift_pos - 3]) == 'T':

                    grna_str += (str(grna) + " ")  # write done the gRNA string
                    grna_pos_str += (str(chr_pos + shift_pos - seq_len) + " ")  # write down PAM sequence position
                    optimal_num += 1
            else:
                grna_str += (str(grna) + " ")
                grna_pos_str += (str(chr_pos + shift_pos - seq_len) + " ")  # write down PAM sequence position

            pam_num += 1  # count the pam num
        shift_pos += 1

    return grna_str, grna_pos_str, pam_num, optimal_num

# sacas9 type for C->G mutation. Pam sequence: NCYYA or NCYYN  Y = T or C
def get_pam_sequence_tc_sacas9(fasta_data, chr_num, chr_pos, shift_pos = 14, seq_len = 20, is_optimal = True):

    pam_num = 0  # count the # of pam sequences
    optimal_num = 0  # count the # of optimal pam sequences
    grna_str = ""  # gRNA string for stored Pam Sequences with space separation
    grna_pos_str = ""     # Pam Sequence Postion for stored Pam Sequence positions with space separation

    pam_checking_times = 4  # the total number of pam searching window shifts
    pam_checking_pos = shift_pos + pam_checking_times  # the max shift position

    while shift_pos < pam_checking_pos:
        isPam = False
        pam_str = fasta_data['chr' + str(chr_num)][chr_pos - shift_pos:chr_pos - shift_pos + 4] # get 4 characters upstream

        # check whether it is PAM NGRRT (NCTTA, NCTCA, NCCTA, NCCCA, NCTTN, NCTCN, NCCTN, NCCCN
        if (str(pam_str) == 'CTTA') or (str(pam_str) == 'CTCA') or (str(pam_str) == 'CCTA') or (str(pam_str) == 'CCCA') \
                or (str(pam_str[0:3]) == 'CTT') or (str(pam_str[0:3]) == 'CTC') or (str(pam_str[0:3]) == 'CCT') or (str(pam_str[0:3]) == 'CCC'):
            isPam = True

        if isPam:
            # Get the PAM sequence from the PAM position upstream 20 characters
            grna = fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - seq_len:chr_pos - shift_pos]

            # check whether it is optimal
            if is_optimal:
                if str(fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - 4:chr_pos - shift_pos - 3]) == 'A' or \
                        str(fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - 4:chr_pos - shift_pos - 3]) == 'T':
                    grna_seq = Seq(str(grna))  # convert it to BioSequence
                    grna_tc = grna_seq.reverse_complement()
                    grna_str += (str(grna_tc) + " ")  # write done the gRNA string
                    grna_pos_str += (str(chr_pos - shift_pos -seq_len) + " ")  # write down PAM sequence position
                    optimal_num += 1
            else:
                grna_seq = Seq(str(grna))
                grna_tc = grna_seq.reverse_complement()
                grna_str += (str(grna_tc) + " ")
                grna_pos_str += (str(chr_pos - shift_pos -seq_len) + " ")  # write down PAM sequence position

            pam_num += 1  # count the pam num
        shift_pos += 1

    return grna_str, grna_pos_str, pam_num, optimal_num


# nmecas9 type for A-G mutation. Pam sequence: NNNNGATT
def get_pam_sequence_ag_nmecas9(fasta_data, chr_num, chr_pos, shift_pos = 14, seq_len = 20, is_optimal = True):  # cas9 is NNNNGATT

    pam_num = 0  # count the # of pam sequences
    optimal_num = 0  # count the # of optimal pam sequences
    grna_str = ""  # gRNA string for stored Pam Sequences with space separation
    grna_pos_str = ""     # Pam Sequence Postion for stored Pam Sequence positions with space separation

    pam_checking_times = 4  # the total number of pam searching window shifts
    pam_checking_pos = shift_pos + pam_checking_times  # the max shift position

    while shift_pos < pam_checking_pos:
        isPam = False
        pam_str = fasta_data['chr' + str(chr_num)][chr_pos + shift_pos:chr_pos + shift_pos + 7]  # get 7 characters downstram

        # check whether it is PAM N-NNNGATT
        if (str(pam_str[3:]) == 'GATT') :
            isPam = True

        if isPam:
            # Get the PAM sequence from the PAM position upstream 20 characters
            grna = fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - seq_len:chr_pos + shift_pos]

            # check whether it is optimal
            if is_optimal:
                if str(fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 4:chr_pos + shift_pos - 3]) == 'A' or str(
                        fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 4:chr_pos + shift_pos - 3]) == 'T':
                    grna_str += (str(grna) + " ")  # write done the gRNA string
                    grna_pos_str += (str(chr_pos + shift_pos -seq_len) + " ")  # write down PAM sequence position
                    optimal_num += 1
            else:
                grna_str += (str(grna) + " ")
                grna_pos_str += (str(chr_pos + shift_pos -seq_len) + " ")  # write down PAM sequence position

            pam_num += 1  # count the pam num
        shift_pos += 1

    return grna_str, grna_pos_str, pam_num, optimal_num

# rmecas9 type for C->G mutation. Pam sequence: N-NNNCTAA
def get_pam_sequence_tc_nmecas9(fasta_data, chr_num, chr_pos, shift_pos = 14, seq_len = 20, is_optimal = True):

    pam_num = 0  # count the # of pam sequences
    optimal_num = 0  # count the # of optimal pam sequences
    grna_str = ""  # gRNA string for stored Pam Sequences with space separation
    grna_pos_str = ""     # Pam Sequence Postion for stored Pam Sequence positions with space separation

    pam_checking_times = 4  # the total number of pam searching window shifts
    pam_checking_pos = shift_pos + pam_checking_times  # the max shift position

    while shift_pos < pam_checking_pos:
        isPam = False
        pam_str = fasta_data['chr' + str(chr_num)][chr_pos - shift_pos:chr_pos - shift_pos + 7]   # get 7 characters upstram

        # check whether it is PAM N-NNNCTAA
        if str(pam_str[3:]) == 'CTAA':
            isPam = True

        if isPam:
            # Get the PAM sequence from the PAM position upstream 20 characters
            grna = fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - seq_len:chr_pos - shift_pos]

            # check whether it is optimal
            if is_optimal:
                if str(fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - 4:chr_pos - shift_pos - 3]) == 'A' or \
                        str(fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - 4:chr_pos - shift_pos - 3]) == 'T':
                    grna_seq = Seq(str(grna))  # convert it to BioSequence
                    grna_tc = grna_seq.reverse_complement()
                    grna_str += (str(grna_tc) + " ")  # write done the gRNA string
                    grna_pos_str += (str(chr_pos - shift_pos -seq_len) + " ")  # write down PAM sequence position
                    optimal_num += 1
            else:
                grna_seq = Seq(str(grna))  # convert it to BioSequence
                grna_tc = grna_seq.reverse_complement()
                grna_str += (str(grna_tc) + " ")
                grna_pos_str += (str(chr_pos - shift_pos -seq_len) + " ")  # write down PAM sequence position

            pam_num += 1  # count the pam num
        shift_pos += 1

    return grna_str, grna_pos_str, pam_num, optimal_num

# cjcas9 type for A-G mutation. Pam sequence: N-NNNRYAC  R = A or G; Y = C or T
# Pam sequence: N-NNNACAC, N-NNNATAC, N-NNNGCAC, N-NNNGTAC
def get_pam_sequence_ag_cjcas9(fasta_data, chr_num, chr_pos, shift_pos = 14, seq_len = 20, is_optimal = True):

    pam_num = 0  # count the # of pam sequences
    optimal_num = 0  # count the # of optimal pam sequences
    grna_str = ""  # gRNA string for stored Pam Sequences with space separation
    grna_pos_str = ""     # Pam Sequence Postion for stored Pam Sequence positions with space separation

    pam_checking_times = 4  # the total number of pam searching window shifts
    pam_checking_pos = shift_pos + pam_checking_times  # the max shift position

    while shift_pos < pam_checking_pos:
        isPam = False
        pam_str = fasta_data['chr' + str(chr_num)][chr_pos + shift_pos:chr_pos + shift_pos + 7]   # get 7 characters downsteam

        # check whether it is PAM N-NNNRYAC (N-NNNACAC, N-NNNATAC, N-NNNGCAC, N-NNNGTAC
        if (str(pam_str[3:]) == 'ACAC') or (str(pam_str[3:]) == 'ATAC') or \
                (str(pam_str[3:]) == 'GCAC') or (str(pam_str[3:]) == 'GTAC') :
            isPam = True

        if isPam:
            # Get the PAM sequence from the PAM position upstream 20 characters
            grna = fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 20:chr_pos + shift_pos]
            # check whether it is optimal
            if is_optimal:
                if str(fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 4:chr_pos + shift_pos - 3]) == 'A' or str(
                        fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 4:chr_pos + shift_pos - 3]) == 'T':
                    grna_str += (str(grna) + " ")  # write done the gRNA string
                    grna_pos_str += (str(chr_pos + shift_pos -seq_len) + " ")  # write down PAM sequence position
                    optimal_num += 1
            else:
                grna_str += (str(grna) + " ")
                grna_pos_str += (str(chr_pos + shift_pos -seq_len) + " ")  # write down PAM sequence position

            pam_num += 1  # count the pam num
        shift_pos += 1

    return grna_str, grna_pos_str, pam_num, optimal_num

# cjcas9 type for T->C mutation. Pam sequence: N-NNNYRTG  R = A or G; Y = C or T
# Pam sequence: N-NNNTGTG, N-NNNTATG, N-NNNCGTG, N-NNNCATG
def get_pam_sequence_tc_cjcas9(fasta_data, chr_num, chr_pos, shift_pos = 14, seq_len = 20, is_optimal = True):

    pam_num = 0  # count the # of pam sequences
    optimal_num = 0  # count the # of optimal pam sequences
    grna_str = ""  # gRNA string for stored Pam Sequences with space separation
    grna_pos_str = ""     # Pam Sequence Postion for stored Pam Sequence positions with space separation

    pam_checking_times = 4  # the total number of pam searching window shifts
    pam_checking_pos = shift_pos + pam_checking_times  # the max shift position

    while shift_pos < pam_checking_pos:
        isPam = False
        pam_str = fasta_data['chr' + str(chr_num)][chr_pos - shift_pos:chr_pos - shift_pos + 7] # get 7 characters upstream

        # check whether it is PAM N-NNNYRAC (N-NNNTGAC, N-NNNTAAC, N-NNNCGAC, N-NNNCAAC
        if (str(pam_str[3:]) == 'TGTG') or (str(pam_str[3:]) == 'TATG') or \
                (str(pam_str[3:]) == 'CGTG') or (str(pam_str[3:]) == 'CATG'):
            isPam = True

        if isPam:
            # Get the PAM sequence from the PAM position upstream 20 characters
            grna = fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - seq_len:chr_pos - shift_pos]

            # check whether it is optimal
            if is_optimal:
                if str(fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - 4:chr_pos - shift_pos - 3]) == 'A' or \
                        str(fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - 4:chr_pos - shift_pos - 3]) == 'T':
                    grna_seq = Seq(str(grna))  # convert it to BioSequence
                    grna_tc = grna_seq.reverse_complement()
                    grna_str += (str(grna_tc) + " ")  # write done the gRNA string
                    grna_pos_str += (str(chr_pos - shift_pos -seq_len) + " ")  # write down PAM sequence position
                    optimal_num += 1
            else:
                grna_seq = Seq(str(grna))  # convert it to BioSequence
                grna_tc = grna_seq.reverse_complement()
                grna_str += (str(grna_tc) + " ")
                grna_pos_str += (str(chr_pos - shift_pos -seq_len) + " ")  # write down PAM sequence position

            pam_num += 1  # count the pam num
        shift_pos += 1

    return grna_str, grna_pos_str, pam_num, optimal_num

# stcas9 type for A-G mutation. Pam sequence: N-NAGAAW  W = A or T;
# Pam sequence: N-NAGAAA, N-NAGAAT
def get_pam_sequence_ag_stcas9(fasta_data, chr_num, chr_pos, shift_pos = 14, seq_len = 20, is_optimal = True):

    pam_num = 0  # count the # of pam sequences
    optimal_num = 0  # count the # of optimal pam sequences
    grna_str = ""  # gRNA string for stored Pam Sequences with space separation
    grna_pos_str = ""     # Pam Sequence Postion for stored Pam Sequence positions with space separation

    pam_checking_times = 4  # the total number of pam searching window shifts
    pam_checking_pos = shift_pos + pam_checking_times  # the max shift position

    while shift_pos < pam_checking_pos:
        isPam = False
        pam_str = fasta_data['chr' + str(chr_num)][chr_pos + shift_pos:chr_pos + shift_pos + 6]   # get 6 characters downsteam

        # check whether it is PAM N-NAGAAW (N-NAGAAA, N-NAGAAT)
        if (str(pam_str[1:]) == 'AGAAA') or (str(pam_str[1:]) == 'AGAAT'):
            isPam = True

        if isPam:
            # Get the PAM sequence from the PAM position upstream 20 characters
            grna = fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - seq_len:chr_pos + shift_pos]

            # check whether it is optimal
            if is_optimal:
                if str(fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 4:chr_pos + shift_pos - 3]) == 'A' or str(
                        fasta_data['chr' + str(chr_num)][chr_pos + shift_pos - 4:chr_pos + shift_pos - 3]) == 'T':
                    grna_str += (str(grna) + " ")  # write done the gRNA string
                    grna_pos_str += (str(chr_pos + shift_pos -seq_len) + " ")  # write down PAM sequence position
                    optimal_num += 1
            else:
                grna_str += (str(grna) + " ")
                grna_pos_str += (str(chr_pos + shift_pos -seq_len) + " ")  # write down PAM sequence position

            pam_num += 1  # count the pam num
        shift_pos += 1

    return grna_str, grna_pos_str, pam_num, optimal_num

# stcas9 type for C->T mutation. Pam sequence: N-NTCTTT, N-NTCTTA
def get_pam_sequence_tc_stcas9(fasta_data, chr_num, chr_pos, shift_pos = 14, seq_len = 20, is_optimal = True):

    pam_num = 0  # count the # of pam sequences
    optimal_num = 0  # count the # of optimal pam sequences
    grna_str = ""  # gRNA string for stored Pam Sequences with space separation
    grna_pos_str = ""     # Pam Sequence Postion for stored Pam Sequence positions with space separation

    pam_checking_times = 4  # the total number of pam searching window shifts
    pam_checking_pos = shift_pos + pam_checking_times  # the max shift position

    while shift_pos < pam_checking_pos:
        isPam = False
        pam_str = fasta_data['chr' + str(chr_num)][chr_pos - shift_pos:chr_pos - shift_pos + 6] # get 6 characters upstream

        # check whether it is N-NTCTTT, N-NTCTTA
        if (str(pam_str[1:]) == 'TCTTT') or (str(pam_str[1:]) == 'TCTTA'):
            isPam = True

        if isPam:
            # Get the PAM sequence from the PAM position upstream 20 characters
            grna = fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - seq_len:chr_pos - shift_pos]

            # check whether it is optimal
            if is_optimal:
                if str(fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - 4:chr_pos - shift_pos - 3]) == 'A' or \
                        str(fasta_data['chr' + str(chr_num)][chr_pos - shift_pos - 4:chr_pos - shift_pos - 3]) == 'T':
                    grna_seq = Seq(str(grna))  # convert it to BioSequence
                    grna_tc = grna_seq.reverse_complement()
                    grna_str += (str(grna_tc) + " ")  # write done the gRNA string
                    grna_pos_str += (str(chr_pos - shift_pos -seq_len) + " ")  # write down PAM sequence position
                    optimal_num += 1
            else:
                grna_seq = Seq(str(grna))  # convert it to BioSequence
                grna_tc = grna_seq.reverse_complement()
                grna_str += (str(grna_tc) + " ")
                grna_pos_str += (str(chr_pos - shift_pos -seq_len) + " ")  # write down PAM sequence position

            pam_num += 1  # count the pam num
        shift_pos += 1

    return grna_str, grna_pos_str, pam_num, optimal_num
