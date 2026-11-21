#!/usr/bin/env bash
# 目的：依次运行公开的 GSE274546 原始矩阵 QC 与独立重注释。
# 输入：沿用 MILLER_IM_DATA_ROOT 环境变量；未设置时由下游使用仓库 data/。
# 输出：沿用 write/36_gse274546_raw_independent_reannotation/ 及下游既有图件目录。
# 运行：bash scripts/03d_GSE274546运行入口.sh。
# 运行前激活已安装依赖的Python环境；python3由当前PATH提供。

set -euo pipefail

# 先在当前脚本所属仓库执行公开 QC，R 入口仍接收仓库根目录。
root="$(cd "$(dirname "$0")/.." && pwd)"
Rscript "$root/scripts/03b_GSE274546质控.R" "$root"

# 重注释消费上一步 QC 输出；其余未公开步骤不属于本发布入口。
python3 "$root/scripts/03c_GSE274546独立重注释.py"
