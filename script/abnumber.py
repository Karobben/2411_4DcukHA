import os
import pandas as pd
from abnumber import Chain
from Bio import SeqIO
from Bio.Seq import Seq
from parallelbar import progress_map

def SeqNumbering(input):
    seq_name, seq_in = input
    try:
        anno_tmp = Chain(str(Seq(seq_in).translate()), scheme='imgt')
    except:
        try:
            anno_tmp = Chain(str(Seq(seq_in[1:]).translate()), scheme='imgt')
        except:
            try:
                anno_tmp = Chain(str(Seq(seq_in[2:]).translate()), scheme='imgt')
            except:
                return None
    anno_result = {
            'seq_id': seq_name,
            'input': seq_in,
            'seq': anno_tmp.seq, 
            'chain_type': anno_tmp.chain_type, 
            'cdr1': anno_tmp.cdr1_seq,
            'cdr2': anno_tmp.cdr2_seq,
            'cdr3': anno_tmp.cdr3_seq,
            'fr1': anno_tmp.fr1_seq,
            'fr2': anno_tmp.fr2_seq,
            'fr3': anno_tmp.fr3_seq,
            'fr4': anno_tmp.fr4_seq,
            'tail': anno_tmp.tail
            }
    return anno_result

# read the sequences from file
def FA2DIC(Seq1):
    Seq_dict = {} 
    for seq_record in SeqIO.parse(Seq1, "fasta"):
        Seq_dict.update({seq_record.id:str(seq_record.seq)})
    return Seq_dict

Seq_lst = os.popen('ls all_VDJ/*/outs/all_contig.fasta').read().split()
    
for Seq1 in Seq_lst:
    SAMPLE = Seq1.split('/')[1]
    Seq_dict = FA2DIC(Seq1)
    input = [(i,ii) for i,ii in zip(list(Seq_dict.keys()), list(Seq_dict.values()))]
    res = progress_map(SeqNumbering, input, process_timeout=1.5, n_cpu=60)
    TB = pd.DataFrame([i for i in res if i != None])
    TB.to_csv(f"result/{SAMPLE}.tsv", sep ='\t')
