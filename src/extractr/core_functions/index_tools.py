def print_kmer_right_read_fragments(kmer, kmer2tf, read_length=310, topk=100, split_springs=True):
  pos = kmer2tf.pos(kmer)
  for p in pos[:topk]:
    if split_springs:
      print(kmer2tf.reads[p:p+read_length].split(b"\n")[0].split("~")[0])
    else:
      print(kmer2tf.reads[p:p+read_length].split(b"\n")[0])

def print_kmer_left_read_fragments(kmer, kmer2tf, topk=100, read_length=310, k=23, split_springs=True):
  pos = kmer2tf.pos(kmer)
  for p in pos[:topk]:
    if split_springs:
      print(kmer2tf.reads[p-read_length:p+k].split(b"\n")[-1].split(b"~")[-1])
    else:
      print(kmer2tf.reads[p-read_length:p+k].split(b"\n")[-1])


def get_kmer_right_read_fragments(kmer, kmer2tf, read_length=310, topk=100, split_springs=True):
  pos = kmer2tf.pos(kmer)
  reads = []
  for p in pos[:topk]:
    if split_springs:
      reads.append(kmer2tf.reads[p:p+read_length].split(b"\n")[0].split(b"~")[0])
    else:
      reads.append(kmer2tf.reads[p:p+read_length].split(b"\n")[0])
  return reads

def print_kmer_left_read_fragments(kmer, kmer2tf, topk=100, read_length=310, k=23, split_springs=True):
  pos = kmer2tf.pos(kmer)
  reads = []
  for p in pos[:topk]:
    if split_springs:
      reads.append(kmer2tf.reads[p-read_length:p+k].split(b"\n")[-1].split(b"~")[-1])
    else:
      reads.append(kmer2tf.reads[p-read_length:p+k].split(b"\n")[-1])
  return reads
