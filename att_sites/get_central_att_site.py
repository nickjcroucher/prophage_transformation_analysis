#! python

import os
import sys
from Bio import SeqIO

att_site_length = 25

rc_list = ["att_ssbB","att_tgt","att_MM1","att_comYC","att_pfbA","att_relA","att_ply"]

with open('att_site_sequences.fa','r') as seq_file:
  for record in SeqIO.parse(seq_file,'fasta'):
    if not record.id.endswith('_alt'):
      seq_length = len(str(record.seq))
      start_point = int(seq_length/2-att_site_length/2)
      end_point = int(seq_length/2+att_site_length/2)
      if record.id in rc_list:
        central_seq = str(record.seq.reverse_complement()[start_point:end_point])
      else:
        central_seq = str(record.seq[start_point:end_point])
      print('>' + record.id + '\n' + central_seq)
