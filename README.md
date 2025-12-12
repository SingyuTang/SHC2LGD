# SHC2LGD

## ITSG-LGD: Along-Orbit Line-of-Sight Gravity Difference Analysis

## 项目简介 (Overview)

本项目实现了基于 **GRACE-FO** 卫星轨道和 **ITSG (Institute of Geodesy at Graz University of Technology)** 每日重力场解算（Daily Solutions）计算 **视线方向重力差 (Line-of-Sight Gravity Difference, LGD)** 的算法。

该方法参考了 *Ghobadi-Far et al. (2022)* 在 JGR: Solid Earth 发表的论文 *"Along-Orbit Analysis of GRACE Follow-On Inter-Satellite Laser Ranging Measurements for Sub-Monthly Surface Mass Variations"*。与论文主要分析 LRI 观测残差不同，本项目侧重于从 **ITSG Daily Gravity Field Models** 中合成高频 LGD 信号，用于分析亚月（Sub-monthly）尺度的质量变化（如洪水、海洋高频信号等）。

## 主要功能 (Features)

1. **高效 LGD 计算 (sh2lgd.py)**:
   - 基于球谐系数 (SHC) 和卫星笛卡尔坐标，计算沿轨 LGD。
   - 采用**向量化 (Vectorized)** 和**递归 (Recursive)** 算法计算勒让德函数及其导数，极大提高了计算大量轨道点时的效率。
   - 支持**扣除静态背景场**（如 GGM/GSM），提取时变信号。
   - 支持区域筛选（只计算经过特定经纬度区域的轨道段）。
2. **信号时频分析与滤波 (MATLAB Scripts)**:
   - 使用连续小波变换 (**CWT**, Morse wavelet) 对 LGD 时间序列进行分析。
   - 提取特定频段（如 1-30 mHz）的信号，去除噪声并重构高频重力信号。
3. **可视化工具 (Python Scripts)**:
   - **单日对比**: 绘制滤波前后的 LGD 空间分布对比。
   - **多日序列**: 绘制特定区域内的 LGD 随时间演变的 "瀑布图" (Waterfall/Offset Plot)，直观展示质量异常的移动和演变。

更多信息查看[博客-项目介绍](https://singyutang.github.io/2025/12/10/%E5%9F%BA%E4%BA%8EITSG-Grace-operational-Kalman-n40%E5%8D%95%E6%97%A5%E7%90%83%E5%8D%8F%E4%BA%A7%E5%93%81%E8%AE%A1%E7%AE%97LGD/)。

## 项目结构 (Project Structure)

codeTextdownloadcontent_copyexpand_less

```
ITSG-LGD-Analysis/
├── sh2lgd.py                        # [核心] Python主程序：读取ITSG模型和轨道，计算LGD并保存为.mat
├── S03cwt_fliter_lgd_ra.m           # [处理] MATLAB脚本：批量调用CWT滤波函数
├── process_cwt_fliter_signal_data.m # [算法] MATLAB函数：执行CWT变换、频带提取和信号重构
├── plot_multi_lgd_cross_over_area.py# [绘图] Python脚本：绘制多日LGD序列在特定区域的变化
├── plot_single_lgd.py               # [绘图] Python脚本：绘制单日滤波前后对比
├── S02compute_grace_lgd.py          # [依赖] 包含OrbitLoader类，用于加载轨道数据 (需确保在路径中)
├── S05plot_lgd_ra_cwt_filter.py     # [依赖] 包含区域筛选工具函数 (需确保在路径中)
└── results/                         # [输出] 存放计算结果 (.mat) 和 图片 (.png)
```

## 环境依赖 (Requirements)

### Python

- Python 3.8+

- NumPy

- SciPy

- Matplotlib

- 自定义依赖 (需确保主要脚本同级目录下存在以下模块):

  ​	-  S02compute_grace_lgd: 用于加载卫星轨道 (OrbitLoader)。

  ​	-  S05plot_lgd_ra_cwt_filter: 用于筛选特定区域的轨道 (filter_complete_tracks_passing_region)。

### MATLAB

- MATLAB R2019b 或更高版本

## 数据准备 (Data Preparation)

在运行代码前，请确保以下数据已下载并按目录结构存放：

1. **ITSG Daily Solutions (.gfc)**:
   - 来源: [ITSG-Grace_operational_Kalman_daily_n40](https://ftp.tugraz.at/outgoing/ITSG/GRACE/ITSG-Grace_operational/daily_kalman/daily_n40/)
   - 存放路径示例: BASE_DIR/itsg_dataset/YYYY-MM/
2. **Static Background Model**:
   - 用于扣除平均场（如GOCO06s或GRACE-FO Level2 GSM产品）。
   - 存放路径示例: `'BASE_DIR\grace_products\GSM-2_2020122-2020152_GRFO_UTCSR_BB01_0603'`
3. **GRACE-FO Orbit Data (GNV1B)**:
   - 包含卫星的位置。代码中通过 OrbitLoader 类加载 GROOPS 工作区数据，你需要根据实际情况修改 `load_orbit` 函数以适配你的轨道文件格式（如 GNV1B等）。
   - 如果按照博客[详细介绍利用GRACE1B多日数据计算LGD工作流程二_基于LRI1B多日数据](https://singyutang.github.io/2025/11/10/详细介绍利用GRACE1B多日数据计算LGD工作流程二-基于LRI1B多日数据/)已经创建了工作根目录`workdir`并进行了相关处理，只需要将本项目文件全部拷贝到该目录下即可运行，因为本项目中的 `load_orbit` 函数默认读取`workdir/gracefo_dataset`路径下的GNV1B轨道文件。

*注：BASE_DIR为GROOPS工作目录，如`BASE_DIR = r'G:\GROOPS\PNAS2020Workspace'`*

## 使用指南 (Usage)

### 第一步：计算合成 LGD (sh2lgd.py)

修改 sh2lgd.py 中的配置区域：

```python
# 配置示例
BASE_DIR = r'G:\GROOPS\PNAS2020Workspace'		# 工作空间根目录
STATIC_FILE_REL = r'grace_products\GSM-2_2020122-2020152_GRFO_UTCSR_BB01_0603'		# 背景场文件 (相对于 BASE_DIR 或绝对路径)
START_DATE = '2020-06-01'	# 处理时间范围起始时间
END_DATE = '2020-08-30'		# 处理时间范围终止时间
MAX_DEGREE = 40  # 根据ITSG模型的阶数设定
REGION_CFG = {
        'lon': (88, 92),  # 感兴趣区域经度范围
        'lat': (22, 26),  # 感兴趣区域纬度范围
        'lat_limit': (-80.0, 80.0),  # 轨道纬度限制
        'dir': 'asc'  # 轨道方向 'asc' (升轨) 或 'desc' (降轨)
    }
```

运行脚本：

```bash
python sh2lgd.py
```

*输出：results/time-lgd-YYYY-MM-DD.mat (包含时间和原始 LGD 数据)*

### 第二步：信号滤波 (S03cwt_fliter_lgd_ra.m)

在 MATLAB 中打开并运行 S03cwt_fliter_lgd_ra.m。

- **配置**：设置 start_datenum, end_datenum 和频率范围（f_high，例如 30 mHz）。
- **功能**：脚本会自动读取上一步生成的 .mat 文件，应用 CWT 滤波保留高频/亚月信号。
- **输出**：results/cwt_time-lgd-YYYY-MM-DD.mat



*注意：这一步必须运行，但是结果不会使用滤波后的数据，因为作者偷懒没有修改结果可视化部分的代码（从其他项目复制过来的代码逻辑），所以包含本部分运行结果的的逻辑部分。*

### 第三步：结果可视化

**1. 多日序列区域分析 (plot_multi_lgd_cross_over_area.py)**

用于展示特定区域（如孟加拉湾、亚马逊流域等）随时间的重力变化。

```python
# 修改日期列表和区域范围
date_list = ['2020-06-04', '2020-06-10', ...]
lon_range = (88, 92)
lat_range = (22, 26)
lat_limit = (-80.0, 80.0)  # 轨道延申纬度范围
orbit_direction = 'asc'  # 轨道方向
```

运行：

```bash
python plot_multi_lgd_cross_over_area.py
```

**2. 单日对比分析 (plot_single_lgd.py)**

查看特定日期的原始信号与滤波后信号的差异。

```python
date_str = '2020-07-07'
lon_range = (88, 92)  # 经度范围（度）
lat_range = (22, 26)  # 纬度范围（度）
lat_limit = (-80.0, 80.0)  # 轨道延申纬度范围
orbit_direction = 'asc'  # 轨道方向
```

运行：

```bash
python plot_single_lgd.py
```

## 算法原理 (Methodology)

### 1. LGD 计算

LGD (Line-of-Sight Gravity Difference) 定义为两颗卫星位置处的引力矢量在视线方向投影的差值。
$$
L G D = ( \overrightarrow { g _ { 2 } } - \overrightarrow { g _ { 1 } } ) \cdot \overrightarrow { e _ { 1 2 } }
$$
其中：

$$ \overrightarrow { g _ { 1 } }$$ 是由 ITSG 球谐系数 $$ ( C _ { n m } , S _ { n m } )$$ 计算得到的重力加速度矢量。本项目使用了**递归算法**直接在轨道高度计算引力梯度。

$$ \overrightarrow { e _ { 1 2 } }$$ 是两颗卫星连线的单位向量。

### 2. ITSG 模型应用

不同于论文中使用 Level-2 月平均数据，本项目利用 ITSG 提供的**每日解 (Daily Solutions)**。通过减去长期静态场，我们能够合成出包含高频瞬变信号（Transient signals）的 LGD 序列。

### 3. CWT 滤波

为了分离特定频段的地球物理信号（去除极低频的轨道误差或极高频的噪声），采用 Morse 小波进行连续小波变换，通过重构特定尺度（对应频率f）的小波系数，提取目标频段信号。

## 参考文献 (References)

1.Ghobadi-Far, Khosro, Shin-Chan Han, Christopher M. McCullough, David N. Wiese, Richard D. Ray, Jeanne Sauber, Linus Shihora, and Henryk Dobslaw. 2022. “Along-Orbit Analysis of GRACE Follow-On Inter-Satellite Laser Ranging Measurements for Sub-Monthly Surface Mass Variations.” *Journal of Geophysical Research: Solid Earth* 127(2):e2021JB022983. doi:[10.1029/2021JB022983](https://doi.org/10.1029/2021JB022983).

2.Mayer-Gürr, Torsten; Behzadpur, Saniya; Ellmer, Matthias; Kvas, Andreas; Klinger, Beate; Strasser, Sebastian; Zehentner, Norbert (2018): ITSG-Grace2018 - Monthly, Daily and Static Gravity Field Solutions from GRACE. GFZ Data Services. doi.org: [10.5880/ICGEM.2018.003](https://doi.org/10.5880/ICGEM.2018.003).
