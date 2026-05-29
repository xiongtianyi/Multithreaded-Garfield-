#!/bin/bash
###############################################################################
# File: run_simulation.sh
# Description: 自动化批量运行 COMSOL、Garfield (THGEM) 和 ROOT 分析。
# Author: YourName
# Updated by ChatGPT
###############################################################################

#----------------------#
#      配置与预设      #
#----------------------#
set -euo pipefail

COUNT=300                       # 外部循环次数
GARFIELD="./thgem"                # Garfield可执行文件
ANALYSIS="./analysis"           # ROOT分析工具
COMSOL_MODEL="comsol/THGEM.mph" # COMSOL模型文件
INPUT_FILE="output.txt"         # 输入参数文件
LOG_FILE="run_simulation.log"   # 日志文件
COMSOL_OUT_PREFIX="comsol/output"

#----------------------#
#    函数与工具函数    #
#----------------------#
function error_exit() {
  echo "[ERROR] $*" >&2
  exit 1
}

function check_executable() {
  [[ -x "$1" ]] || error_exit "可执行文件 $1 不存在或不可执行。"
}

function check_file_exists() {
  [[ -f "$1" ]] || error_exit "文件 $1 不存在。"
}

#----------------------#
#     主脚本开始       #
#----------------------#
exec > >(tee -a "$LOG_FILE") 2>&1
echo "============ 开始批量模拟脚本 ============"
date "+当前时间：%Y-%m-%d %H:%M:%S"

# 检查依赖项
check_executable "$GARFIELD"
check_executable "$ANALYSIS"
check_file_exists "$COMSOL_MODEL"
check_file_exists "$INPUT_FILE"

# # 调试输出：验证参数数量与值
# echo "[DEBUG] 参数名数量: ${#PARAM_NAMES[@]}"
# echo "[DEBUG] 参数名: ${PARAM_NAMES[@]}"
# echo "[DEBUG] 参数值数量: ${#PARAM_VALUES[@]}"
# echo "[DEBUG] 参数值: ${PARAM_VALUES[@]}"

#----------------------#
#    主循环 (1..COUNT) #
#----------------------#
for (( i = 1; i <= COUNT; i++ )); do
  echo "---------------------------------------------------"
  echo "[INFO] 外部循环 #$i / $COUNT"
  echo "---------------------------------------------------"

  #----------------------#
  # 修正：读取参数名和参数值
  #----------------------#
  # 读取参数名（第一行，逗号分隔）
  IFS=',' read -r -a PARAM_NAMES < <(head -n 1 "$INPUT_FILE")

  # 读取参数值（第二行，空格分隔）
  IFS=' ' read -r -a PARAM_VALUES < <(tail -n 1 "$INPUT_FILE")

  # 检查参数数量是否匹配
  if [[ ${#PARAM_NAMES[@]} -ne ${#PARAM_VALUES[@]} ]]; then
    error_exit "参数名与参数值数量不匹配。"
  fi

  #----------------------#
  #   (A) 运行 COMSOL    #
  #----------------------#
  # 拼接参数名和值（COMSOL格式：逗号分隔）
  param_names_str=$(IFS=','; echo "${PARAM_NAMES[*]}")
  param_values_str=$(IFS=','; echo "${PARAM_VALUES[*]}")

  # 执行 COMSOL
  echo "[INFO] 运行 COMSOL，参数名: $param_names_str"
  echo "[INFO] 运行 COMSOL，参数值: $param_values_str"
  comsol batch \
    -inputfile "$COMSOL_MODEL" \
    -pname "$param_names_str" \
    -plist "$param_values_str" \
    -outputfile "${COMSOL_OUT_PREFIX}_${i}.mph" || {
      error_exit "COMSOL 在外部循环 $i 时运行失败。"
    }

  # 清理临时文件
  rm -f "${COMSOL_OUT_PREFIX}_${i}.mph.recovery" \
        "${COMSOL_OUT_PREFIX}_${i}.mph.status" \
        "${COMSOL_OUT_PREFIX}_${i}.mph"

  #----------------------#
  # (B) 运行 Garfield    #
  #----------------------#
  echo "[INFO] 执行 Garfield，参数：$i"
  "$GARFIELD" "$i" || error_exit "Garfield 运行失败。"

  #----------------------#
  # (C) 运行 ROOT 分析   #
  #----------------------#
  echo "[INFO] 运行 ROOT 分析：$ANALYSIS $i"
  "$ANALYSIS" "$i" || error_exit "ROOT 分析失败。"

  echo "[INFO] 外部循环 #$i 已完成。"
done

echo "=========================================="
echo "[INFO] 脚本执行完毕，所有循环均已完成。"
date "+当前时间：%Y-%m-%d %H:%M:%S"
echo "=========================================="