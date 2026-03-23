#!/bin/bash
# Usage: bash scripts/run_jcvi.sh
set -e

CPU=64
WORKDIR=~/DATA/Fritillaria_cirrhosa/Genome/collinear_jcvi
cd $WORKDIR

#----------------------------------------------------------
# 第一阶段：格式转换与去冗余
#----------------------------------------------------------
echo "Step 1: Data Formatting..."
species=("Fcr" "Lre" "Lsa")

for sp in "${species[@]}"; do
    mkdir -p work/01.format/$sp
    # 1. 提取BED (JCVI风格)
    # 注意：需根据GFF实际内容调整 --type 和 --id_attr
    python3 -m jcvi.formats.gff bed --type=mRNA --key=ID data/${sp}.gff -o work/01.format/${sp}.bed
    
    # 2. 仅保留前12条大染色体 (针对川贝母和百合)
    head -n 12 work/01.format/${sp}.bed > work/01.format/${sp}.ids.txt
    
    # 3. 提取最长蛋白序列 (假设您已有原始pep.fa，或使用gffread提取)
    # 这里以您已有的pep.fa为准
    python3 scripts/prep_data.py longest data/${sp}.pep.fa work/01.format/${sp}.pep
    
    # 4. 生成LEN文件
    python3 scripts/prep_data.py len data/${sp}.gff work/01.format/${sp}.len
done

#----------------------------------------------------------
# 第二阶段：Diamond 比对 (All-against-all)
#----------------------------------------------------------
echo "Step 2: Diamond Alignment..."
cd work/02.blast
# 建立索引
for sp in "${species[@]}"; do
    diamond makedb --in ../01.format/${sp}.pep -d ${sp}
done

# 两两比对 (Fcr vs Lre, Fcr vs Lsa, Lre vs Lsa)
pairs=("Fcr.Lre" "Fcr.Lsa" "Lre.Lsa")
for pair in "${pairs[@]}"; do
    q=${pair%.*}
    s=${pair#*.}
    diamond blastp -d $s.dmnd -q ../01.format/$q.pep -o ${q}.${s}.blast \
        --evalue 1e-5 --max-target-seqs 20 --outfmt 6 --threads $CPU
done

#----------------------------------------------------------
# 第三阶段：JCVI 共线性分析
#----------------------------------------------------------
echo "Step 3: Synteny Detection..."
cd ../03.jcvi
# 链接必要文件到当前目录 (JCVI要求bed和blast同名且在同目录)
for sp in "${species[@]}"; do
    ln -sf ../01.format/${sp}.bed .
    ln -sf ../01.format/${sp}.pep .
done
for pair in "${pairs[@]}"; do
    ln -sf ../02.blast/${pair}.blast .
done

# 运行共线性搜索
python3 -m jcvi.compara.synteny screen Fcr.Lre.blast Fcr.Lre.bed
python3 -m jcvi.compara.synteny screen Fcr.Lsa.blast Fcr.Lsa.bed

# 绘制点图 (Dotplot)
python3 -m jcvi.graphics.dotplot Fcr.Lre.anchors
