## 川贝母与百合属共线性分析流程

针对**巨大基因组**（30Gb - 50Gb）设计的共线性分析工作流，重点解决了 32 位整数溢出、磁盘随机 I/O 瓶颈以及长耗时任务的断点续传问题。

### 1. 💻 计算资源与配置要求
| 资源类型 | 配置要求 | 说明 |
| :--- | :--- | :--- |
| **操作系统** | Ubuntu 22.04 LTS | 核心生产环境 |
| **内存 (RAM)** | ≥ 128 GB (推荐 1T) | 满足 Diamond 大 block-size 比对 |
| **存储空间** | ≥ 500 GB 可用 | 存放索引、比对中间件及 .pep 文件 |
| **关键依赖** | Python 3.8+, Biopython, **pyfaidx**, Diamond, JCVI | `pip install pyfaidx more_itertools` |

---

### 2. 输入数据结构 (Data Preparation)
原始数据通过软链接方式统一存放于 `data/` 目录，并预先生成 `.fai` 索引以加速提取。

* **Fcr (川贝母)**: `Fcr.genome.fa`, `Fcr.gff`
* **Lre (王百合)**: `Lre.genome.fa`, `Lre.gff`
* **Lsa (岷江百合)**: `Lsa.genome.fa`, `Lsa.gff`

---

### 3. 核心用法 (Scripts Usage)

#### **A. 数据预处理脚本 (`prep_data.py`)**
* **功能**：基于 `.fai` 索引的 64 位寻址序列提取、去冗余及三列式长度统计。
* **用法**：
    ```bash
    # 提取蛋白
    python3 scripts/prep_data.py extract <gff> <fasta> <out_pep>
    # 筛选最长转录本
    python3 scripts/prep_data.py longest <in_pep> <out_pep>
    # 生成三列式 LEN 文件 (Name, Length, GeneCount)
    python3 scripts/prep_data.py len <gff> <out_len>
    ```

#### **B. 主控流水线 (`check_and_run.sh`)**
* **功能**：自动化执行格式转换、Diamond 数据库构建、两两比对。内置 **Checkpoint** 机制，存在结果文件时自动跳过。
* **性能优化**：`--block-size 20.0` 利用大内存加速；`--threads 80` 并行计算。
* **启动命令**：
    ```bash
    bash scripts/check_and_run.sh
    ```

---

### 4. 分析逻辑 (Analysis Pipeline)



1.  **Stage 1: Format Conversion**
    * 将 GFF 转换为标准 BED6 格式。
    * 通过 `pyfaidx` 绕过磁盘读取瓶颈，提取 CDS 并翻译为 PEP。
2.  **Stage 2: Homology Search**
    * 使用 Diamond 对蛋白质序列进行 All-against-all 比对。
    * 输出格式为 `outfmt 6` 的比对矩阵。
3.  **Stage 3: Synteny Detection**
    * 利用 JCVI 的 `mcscan` 算法搜索共线性锚点（Anchors）。

---

### 5. 输出文件列表 (Final Outputs)

#### **中间产物 (`work/`)**
* **PEP 文件**：经过去冗余的物种蛋白库（~50MB/物种）。
* **BED 文件**：物种基因坐标参考。
* **LEN 文件**：包含染色体长度与基因密度的三列统计表。
* **BLAST 文件**：物种间同源对原始记录（如 `Fcr.Lre.blast`）。

#### **共线性核心 (`work/03.jcvi/`)**
* **`.anchors`**：核心锚点文件，记录了共线性块内的基因对。
* **`.simple`**：共线性块的简化摘要，用于绘制布局图。
* **`.pdf`**：
    * `Fcr.Lre.pdf`: 散点图 (Dotplot)，展示全基因组共线性趋势。
    * `layout.pdf`: 线性布局图，展示染色体间的物理连接关系。

---
