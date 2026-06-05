# THGEM Charging-up 多线程模拟

本目录用于 THGEM 充电效应模拟。核心流程是：COMSOL 计算电场，`gem` 用 Garfield++ 做多线程雪崩和离子漂移，`analysis` 根据电子/离子端点更新介质表面的累计电荷，下一轮 COMSOL 再读入累计电荷重新计算电场。

## 目标环境

推荐在以下环境中运行：

| 项目 | 版本/要求 |
|------|-----------|
| 操作系统 | Ubuntu 20.04 LTS |
| 编译器 | GCC/G++ 9，按 C++11 编译 |
| CMake | 3.12 或更高 |
| ROOT | 6.22.06 |
| Garfield++ | 修改了Random.cc/hh 降低了多线程瓶颈 |
| COMSOL | COMSOL Multiphysics 6.3 |

`CMakeLists.txt` 已显式设置：

```cmake
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
```

## 目录结构

```text
Charging-up/
├── gem.C                # 多线程 Garfield++ 主模拟程序
├── analysis.c           # 从 ROOT 端点计算并累积表面电荷
├── output.c             # 可选：增益随迭代轮数变化的后处理
├── position.c           # 可选：电子/离子端点三维可视化
├── simulate.sh          # COMSOL -> gem -> analysis 自动迭代脚本
├── CMakeLists.txt       # 构建配置
├── output.txt           # COMSOL 读取的累计表面电荷参数
├── comsol/
│   └── THGEM.mph        # COMSOL 6.3 模型文件
└── simdata/
    ├── THGEM.mphtxt     # COMSOL 导出的网格
    ├── THGEM.txt        # COMSOL 导出的电势数据
    └── dielectrics.dat  # 材料域定义
```

构建后会在 `build/` 下生成 `result/`、`fig/`、`simdata/` 和 `comsol/`。

## 环境变量

先加载 ROOT 6.22.06：

```bash
source /path/to/root-6.22.06/bin/thisroot.sh
root-config --version
```

再加载 Garfield++：

```bash
source /path/to/garfieldpp/install/share/Garfield/setupGarfield.sh
echo "$GARFIELD_INSTALL"
```

如果 Garfield++ 还没有安装，可在仓库根目录执行：

```bash
cd /path/to/garfieldpp
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/garfieldpp/install ..
make -j"$(nproc)"
make install
source /path/to/garfieldpp/install/share/Garfield/setupGarfield.sh
```

确保 COMSOL 6.3 命令行可用：

```bash
export PATH=/usr/local/comsol63/multiphysics/bin:$PATH
comsol -version
```

建议把 ROOT、Garfield++ 和 COMSOL 的环境加载命令写入运行机器的 `~/.bashrc`，避免批量迭代时找不到库或 `comsol` 命令。

## 编译

进入本目录后创建独立构建目录：

```bash
cd /path/to/garfieldpp/Examples/Charging-up
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j"$(nproc)"
```

编译成功后应得到：

```text
build/gem
build/analysis
build/output
build/position
```

同时 `simdata/`、`comsol/`、`output.txt` 和 `simulate.sh` 会被复制到 `build/`。

## gem 主程序

`gem` 是基础多线程模拟程序，负责：

- 读取 `simdata/THGEM.mphtxt`、`simdata/THGEM.txt` 和 `simdata/dielectrics.dat`
- 设置 Ne/CH4 = 95/5 气体，温度 293.15 K，压力 760 Torr
- 使用 `AvalancheMicroscopic` 计算电子雪崩
- 使用 `AvalancheMC` 漂移离子
- 将电子/离子端点写入 `result/resultN.root`
- 使用原子事件计数器把事件分配给多个线程

运行格式：

```bash
./gem [N] [chunk_size] [num_threads] [total_events]
```

| 参数 | 含义 | 默认值 |
|------|------|--------|
| `N` | 迭代轮数，输出 `result/resultN.root` | 未给出时输出 `result/result.root` |
| `chunk_size` | 每次从事件队列领取的事件数 | 1 |
| `num_threads` | 工作线程数 | 32 |
| `total_events` | 总事件数 | 10000 |

示例：

```bash
cd build
./gem 1 1 32 10000
```

## ROOT 输出结构

`gem` 只写一个名为 `Tree` 的 TTree：

| 分支 | 类型 | 每组数据含义 |
|------|------|--------------|
| `e1hit` | `vector<double>` | 电子起点：`x, y, z, t, E` |
| `e2hit` | `vector<double>` | 电子终点：`x, y, z, t, E` |
| `i1hit` | `vector<double>` | 离子起点：`x, y, z, t` |
| `i2hit` | `vector<double>` | 离子漂移终点：`x, y, z, t` |
| `evt` | `int` | 事件编号 |
| `ne` | `int` | 雪崩电子数 |
| `ni` | `int` | 雪崩离子数 |

这些分支由 `analysis.c` 读取，用于统计介质厚度方向上的电子/离子终点数。

## 电荷累积逻辑

`analysis` 读取 `result/resultN.root` 后，会筛选 THGEM 介质板范围内的端点：

```text
z in [induce + metal, induce + metal + thick]
```

当前几何参数为：

| 参数 | 值 | 单位 |
|------|----|------|
| `pitch` | 0.1 | cm |
| `metal` | 0.0012 | cm |
| `thick` | 0.08 | cm |
| `induce` | 0.2 | cm |
| `driftz` | 0.5 | cm |

`analysis` 将介质厚度方向分为 22 个 bin，并把每个 bin 的净电荷累加到 `output.txt`：

```text
delta_charge = (ion_count - electron_count) * 30
new_charge = old_charge + delta_charge
```

`output.txt` 的格式为：

```text
n1,n2,n3,...,n22
0 0 0 ... 0
```

第一行是 COMSOL 参数名，第二行是每个参数对应的累计电荷值。COMSOL 模型需要用这些参数更新表面电荷边界条件。

## 手动运行一轮

```bash
cd build

# 第 1 轮 Garfield++ 模拟
./gem 1 1 32 10000

# 根据第 1 轮端点更新 output.txt
./analysis 1
```

如果要继续第 2 轮，需要先让 COMSOL 6.3 读取新的 `output.txt` 并重新导出 `simdata/THGEM.txt`、`simdata/THGEM.mphtxt`，然后再运行：

```bash
./gem 2 5 32 10000
./analysis 2
```

## 自动迭代运行

`simulate.sh` 会自动执行：

```text
读取 output.txt
-> comsol batch 更新电场
-> ./gem N
-> ./analysis N 更新 output.txt
-> 进入下一轮
```

运行：

```bash
cd build
chmod +x simulate.sh
./simulate.sh
```

脚本中常用配置：

```bash
COUNT=100
GARFIELD="./gem"
ANALYSIS="./analysis"
COMSOL_MODEL="comsol/THGEM.mph"
INPUT_FILE="output.txt"
LOG_FILE="run_simulation.log"
```

如果需要改变每轮事件数、线程数或 chunk 大小，可把脚本中的 Garfield 调用改成：

```bash
"$GARFIELD" "$i" 1 32 10000
```

## COMSOL 6.3 注意事项

`simulate.sh` 使用如下命令调用 COMSOL：

```bash
comsol batch \
  -inputfile comsol/THGEM.mph \
  -pname "$param_names_str" \
  -plist "$param_values_str" \
  -outputfile "comsol/output_${i}.mph"
```

请确认 `THGEM.mph` 中的参数名与 `output.txt` 第一行完全一致，例如 `n1` 到 `n22`。还需要确认模型中的导出节点会把新的场图写到 `simdata/THGEM.txt` 和 `simdata/THGEM.mphtxt`，否则 `thgem` 会继续读取旧电场。

## 常见检查

检查 ROOT：

```bash
root-config --version
root -b -q
```

检查 Garfield++：

```bash
echo "$GARFIELD_INSTALL"
echo "$CMAKE_PREFIX_PATH"
echo "$LD_LIBRARY_PATH"
```

检查 COMSOL：

```bash
which comsol
comsol -version
```

检查结果文件：

```bash
root -l result/result1.root
Tree->Print()
```

## 输出文件

| 文件/目录 | 内容 |
|-----------|------|
| `result/resultN.root` | 第 N 轮的端点数据，包含 `Tree` |
| `output.txt` | 累计表面电荷，供 COMSOL 下一轮读取 |
| `fig/ZDistribution_N.pdf` | `analysis` 生成的 Z 分布图 |
| `fig/ZDistribution_N.svg` | `analysis` 生成的 Z 分布图 |
| `run_simulation.log` | `simulate.sh` 批量运行日志 |
| `fig/result1.*` | 可选：`output` 生成的增益演化图 |
| `fig/EndPoints3D-N.pdf` | 可选：`position` 生成的端点三维图 |

## 故障排查

如果 `cmake ..` 找不到 ROOT，先确认已经执行：

```bash
source /path/to/root-6.22.06/bin/thisroot.sh
```

如果找不到 Garfield++，先确认已经执行：

```bash
source /path/to/garfieldpp/install/share/Garfield/setupGarfield.sh
```

如果运行时提示找不到共享库，检查：

```bash
echo "$LD_LIBRARY_PATH"
```

如果 `analysis` 找不到 `Tree`，说明 `gem` 没有成功写出 `result/resultN.root`，或运行目录不是 `build/`。

如果每轮电荷没有变化，优先检查 COMSOL 是否真正重新导出了 `simdata/THGEM.txt`。可以用文件时间戳确认：

```bash
ls -lh --time-style=long-iso simdata/THGEM.txt simdata/THGEM.mphtxt
```
