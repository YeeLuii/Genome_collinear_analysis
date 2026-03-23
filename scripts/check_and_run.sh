#!/bin/bash
# filename: scripts/check_and_run.sh

set -e

# --- 1. 工具链检查函数 ---
check_tool() {
    if command -v $1 >/dev/null 2>&1; then
        echo -e "[OK] $1 已安装: $(which $1)"
    else
        echo -e "[ERROR] 未找到 $1。请激活 conda 环境或安装该工具。"
        exit 1
    fi
}

echo "--- 开始环境检查 ---"
tools=("diamond" "python3" "gffread" "bedtools")
for tool in "${tools[@]}"; do
    check_tool $tool
done

# 检查 Python 模块
python3 -c "import jcvi; import Bio" 2>/dev/null || { echo "[ERROR] 缺少 python 模块: jcvi 或 biopython"; exit 1; }
echo "--- 环境检查完毕，开始处理 ---"

# --- 2. 路径设置 ---
WORKDIR=~/DATA/Fritillaria_cirrhosa/Genome/collinear_jcvi
CPU=40 # 针对 112 线程机器的保守配置
cd $WORKDIR

mkdir -p work/01.format work/02.blast work/03.jcvi

# --- 3. 第一阶段：数据格式化 (不再分割染色体) ---
species=("Fcr" "Lre" "Lsa")

for sp in "${species[@]}"; do
    echo ">>> 处理物种: $sp"
    
    # 提取 PEP
    if [ ! -f "work/01.format/${sp}.pep" ]; then
        echo "   [Run] 提取蛋白序列..."
        python3 scripts/prep_data.py extract data/${sp}.gff data/${sp}.genome.fa work/01.format/${sp}.pep.full
        python3 scripts/prep_data.py longest work/01.format/${sp}.pep.full work/01.format/${sp}.pep
        rm -f work/01.format/${sp}.pep.full
    else
        echo "   [Skip] PEP 文件已存在"
    fi
    
    # 生成 BED
    if [ ! -f "work/01.format/${sp}.bed" ]; then
        echo "   [Run] 生成 BED 文件..."
        python3 -m jcvi.formats.gff bed --type=mRNA --key=ID data/${sp}.gff -o work/01.format/${sp}.bed
    fi

    # 生成三列 LEN (始终更新，确保数据最新)
    python3 scripts/prep_data.py len data/${sp}.gff work/01.format/${sp}.len
done

# --- 4. 第二阶段：Diamond 比对 ---
# Diamond 完美支持大基因组索引
cd work/02.blast
for sp in "${species[@]}"; do
    if [ ! -f "${sp}.dmnd" ]; then
        echo ">>> 构建 Diamond 数据库: $sp"
        diamond makedb --in ../01.format/${sp}.pep -d ${sp} --threads $CPU
    fi
done

pairs=("Fcr.Lre" "Fcr.Lsa")
for pair in "${pairs[@]}"; do
    q=${pair%.*}
    s=${pair#*.}
    if [ ! -f "${q}.${s}.blast" ]; then
        echo ">>> 执行蛋白质比对: $q vs $s"
        # 针对 1T 内存增加 block-size 到 20.0 以大幅提升速度
        diamond blastp -d $s.dmnd -q ../01.format/$q.pep -o ${q}.${s}.blast \
            --evalue 1e-10 --max-target-seqs 20 --block-size 20.0 --index-chunks 1 --threads $CPU --outfmt 6
    fi
done

echo "数据预处理与比对已完成。接下来您可以运行 JCVI 绘图。"
echo "分析完成。结果位于 work/02.blast/ "
