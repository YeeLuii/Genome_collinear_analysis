import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

# 尝试导入更高效的 pyfaidx
try:
    from pyfaidx import Fasta
    HAS_PYFAIDX = True
except ImportError:
    HAS_PYFAIDX = False

def extract_pep_from_gff(gff_path, fasta_path, out_path):
    print(f">>> 正在加载基因组索引: {fasta_path} ...")
    if HAS_PYFAIDX:
        # 真正的 .fai 模式，瞬间加载
        genome = Fasta(fasta_path)
    else:
        # 备选方案：Biopython index
        genome = SeqIO.index(fasta_path, "fasta")
    
    gene_struct = defaultdict(list)
    chrom_map = {}
    
    print(f">>> 解析 GFF 文件...")
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'CDS': continue
            
            chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attrs = dict(item.split('=', 1) for item in parts[8].split(';') if '=' in item)
            parent = attrs.get('Parent', '').strip().split(',')[0] # 处理多Parent情况
            
            if parent:
                gene_struct[parent].append((start, end, strand))
                chrom_map[parent] = chrom

    print(f">>> 开始提取序列 (Total mRNA: {len(gene_struct)}) ...")
    count = 0
    with open(out_path, 'w') as out_f:
        for mrna_id, segments in gene_struct.items():
            chrom = chrom_map[mrna_id]
            if chrom not in genome: continue
            
            strand = segments[0][2]
            # 排序逻辑：正链坐标从小到大，负链从大到小
            segments.sort(key=lambda x: x[0], reverse=(strand == '-'))
            
            try:
                # 拼接 CDS
                if HAS_PYFAIDX:
                    # pyfaidx 使用 0-based, end-exclusive
                    full_cds = "".join([genome[chrom][s-1:e].seq for s, e, _ in segments])
                else:
                    full_cds = "".join([str(genome[chrom].seq[s-1:e]) for s, e, _ in segments])
                
                cds_seq = Seq(full_cds)
                if strand == '-':
                    cds_seq = cds_seq.reverse_complement()
                
                pep_seq = cds_seq.translate(to_stop=True)
                if len(pep_seq) > 0:
                    out_f.write(f">{mrna_id}\n{pep_seq}\n")
            except Exception as e:
                pass # 忽略极少数异常序列
            
            count += 1
            if count % 5000 == 0:
                print(f"   已完成 {count} 条转录本...")

def get_longest_peptide(input_pep, output_pep):
    seqs = defaultdict(list)
    for rec in SeqIO.parse(input_pep, "fasta"):
        gene_id = rec.id.rsplit('.', 1)[0]
        seqs[gene_id].append(rec)
    
    results = []
    for gene_id, mrna_list in seqs.items():
        mrna_list.sort(key=lambda x: len(x.seq), reverse=True)
        results.append(mrna_list[0])
    SeqIO.write(results, output_pep, "fasta")

def make_len_file(gff_path, out_path):
    chrom_lens = {}
    chrom_counts = {}
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.split('\t')
            if len(parts) < 5: continue
            chrom, end_pos = parts[0], int(parts[4])
            if chrom not in chrom_lens or end_pos > chrom_lens[chrom]:
                chrom_lens[chrom] = end_pos
            if parts[2] == 'mRNA':
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
    with open(out_path, 'w') as out_f:
        for chrom in sorted(chrom_lens.keys()):
            out_f.write(f"{chrom}\t{chrom_lens[chrom]}\t{chrom_counts.get(chrom, 0)}\n")

if __name__ == "__main__":
    action = sys.argv[1]
    if action == "extract": extract_pep_from_gff(sys.argv[2], sys.argv[3], sys.argv[4])
    elif action == "longest": get_longest_peptide(sys.argv[2], sys.argv[3])
    elif action == "len": make_len_file(sys.argv[2], sys.argv[3])