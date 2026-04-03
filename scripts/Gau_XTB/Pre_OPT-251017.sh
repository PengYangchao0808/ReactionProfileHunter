#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# 0) 脚本目录
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# 1) 可配置参数
NPROCS="${NPROCS:-1}"
GAUSS_BIN="${GAUSS_BIN:-g16}"

# 三类路线（可用环境变量覆盖）
ROUTE_INT_SP="${ROUTE_INT_SP:-#p opt=(calcall,noeigen,nomicro) external='./xtb.sh' nosymm}"
ROUTE_TS="${ROUTE_TS:-#p opt=(calcall,TS,noeigen,nomicro) external='./xtb.sh' nosymm}"
ROUTE_IRC="${ROUTE_IRC:-#p IRC(calcfc,stepsize=5) external='./xtb.sh' nosymm}"

# 未匹配名的回退（默认与 INT/S/P 相同）
ROUTE_FALLBACK="${ROUTE_LINE:-$ROUTE_INT_SP}"

# 2) xtb 源目录与文件清单
XTBSRC_DIR="${XTBSRC_DIR:-$SCRIPT_DIR/xtb}"
XTBF=("extderi" "extderi.f90" "genxyz" "genxyz.f90" "xtb.sh")

# 3) 是否删除脚本目录中的原始 xtb 文件（危险：默认关闭）
DELETE_XTB_FROM_SCRIPT_DIR="${DELETE_XTB_FROM_SCRIPT_DIR:-0}"

# 4) 收集 .gjf
gjfs=( *.gjf )
if (( ${#gjfs[@]} == 0 )); then
  echo "没有 .gjf 文件，退出。"
  exit 0
fi

echo "将按"文件名"自动选择关键词并运行 Gaussian： ${gjfs[*]}"

# 4.1 统一换行符为 LF（避免空行匹配失败）
for f in "${gjfs[@]}"; do
  sed -i 's/\r//' "$f"
done

# 函数：根据文件名选择路线
select_route_from_name() {
  local name="$1"
  # 转大写，便于无关大小写匹配（Bash 4+）
  local u="${name^^}"
  
  # 1) IRC 任务：TS…-IRC 结尾 或 直接以 IRC 开头
  if [[ "$u" =~ ^TS([_-]|[0-9]).*IRC$ || "$u" =~ ^IRC([_-]|[0-9]).* ]]; then
    printf "%s" "$ROUTE_IRC"
  # 2) TS 任务：以 TS 开头，后接连字符/下划线/数字
  elif [[ "$u" =~ ^TS([_-]|[0-9]).* ]]; then
    printf "%s" "$ROUTE_TS"
  # 3) INT/S/P 任务：以 INT/S/P 开头，后接连字符/下划线/数字
  elif [[ "$u" =~ ^(INT|S|P)([-_]|[0-9]).* ]]; then
    printf "%s" "$ROUTE_INT_SP"
  # 4) 其他：回退（默认等同 INT/S/P，可通过环境变量 ROUTE_LINE 覆盖）
  else
    printf "%s" "$ROUTE_FALLBACK"
  fi
}

# 5) 就地重写头部（覆盖至首个空行，保留其后的内容）
for f in "${gjfs[@]}"; do
  # 取不含扩展名的文件名作为"Title"
  base="$(basename "$f" .gjf)"
  route="$(select_route_from_name "$base")"
  
  if ! grep -q -m1 '^$' "$f"; then
    echo "警告：$f 不含 Route 与 Title 之间的空行（格式异常）。将直接覆盖头部并原样拼接。"
  fi
  
  echo "文件: $f  名称: $base  路线: $route"
  tmp="$(mktemp)"
  {
    printf "%%nprocshared=%s\n" "$NPROCS"
    printf "%s\n\n" "$route"
    if grep -q -m1 '^$' "$f"; then
      # 拼接第一个空行之后的内容（不包含该空行）
      sed -n '1,/^$/d; p' "$f"
    else
      # 回退：保留原文件全部内容（可能重复旧头部，建议先修正输入格式）
      cat "$f"
    fi
  } > "$tmp"
  mv -f "$tmp" "$f"
done

# 6) 复制 xtb 所需文件到当前目录
echo "从 $XTBSRC_DIR 复制 xtb 工具到 $PWD"
for x in "${XTBF[@]}"; do
  src="$XTBSRC_DIR/$x"
  if [[ -f "$src" ]]; then
    cp -f "$src" "./$x"
  else
    echo "警告：缺少 $src"
  fi
done
chmod +x ./xtb.sh 2>/dev/null || true

# 7) 逐一运行 Gaussian
for f in "${gjfs[@]}"; do
  out="${f%.gjf}.out"
  echo "运行 $GAUSS_BIN: $f -> $out"
  time "$GAUSS_BIN" < "$f" > "$out"
done

# 8) 可选删除脚本目录原始 xtb 文件（不建议）
if (( DELETE_XTB_FROM_SCRIPT_DIR == 1 )); then
  echo "危险操作：删除脚本目录中的 xtb 原始文件"
  rm -f "${XTBF[@]/#/$SCRIPT_DIR/}" 2>/dev/null || true
  rm -f "${XTBF[@]/#/$SCRIPT_DIR/xtb/}" 2>/dev/null || true
else
  echo "默认不删除脚本目录原件。若必须删除，设置 DELETE_XTB_FROM_SCRIPT_DIR=1 后再运行。"
fi

# 9) 检查 .out 是否正常结束并归档
mkdir -p Finished
for f in "${gjfs[@]}"; do
  out="${f%.gjf}.out"
  if [[ -f "$out" ]] && grep -q "Normal termination of Gaussian" "$out"; then
    mv -f "$f" "$out" Finished/
    echo "完成：$f 正常结束，已移入 Finished/"
  else
    echo "未完成：$f 未正常结束，保留在当前目录以便修改计算条件。"
  fi
done
