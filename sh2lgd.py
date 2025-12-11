import os
import re
import numpy as np
from datetime import datetime, timedelta
from scipy.io import savemat

from S02compute_grace_lgd import OrbitLoader
from S05plot_lgd_ra_cwt_filter import filter_complete_tracks_passing_region


def datetime_to_datenum(d):
    """Python datetime 转 MATLAB datenum"""
    mdn = datetime.toordinal(d) + (d - datetime(d.year, d.month, d.day)).seconds / 86400.0 + 366
    return mdn


class LGDCalculator:
    def __init__(self, GM, R_ref, max_degree):
        """
        初始化计算器
        GM: 地球引力常数 (m^3/s^2)
        R_ref: 参考半径 (m)
        max_degree: 最大截断阶数 L
        """
        self.GM = GM
        self.R_ref = R_ref
        self.L = max_degree

        # --- 优化 1: 预计算 b 系数 ---
        # 维度 [L+2, L+2]，方便直接索引
        self.b1 = np.zeros((self.L + 2, self.L + 2))
        self.b2 = np.zeros((self.L + 2, self.L + 2))
        self.b3 = np.zeros((self.L + 2, self.L + 2))
        self.b4 = np.zeros((self.L + 2, self.L + 2))

        self._precompute_b_coeffs()

    def _precompute_b_coeffs(self):
        """预先计算所有 l, m 的 b 系数"""
        for l in range(2, self.L + 1):
            for m in range(l + 1):
                # b1
                self.b1[l, m] = np.sqrt((l + 1) * (l + 2) * (2 * l + 1) / (2 * (2 * l + 3)))
                # b2
                self.b2[l, m] = np.sqrt((l + m + 1) * (l + m + 2) * (2 * l + 1) / (2 * l + 3))
                # b3 (delta_1m check)
                delta_1m = 1.0 if m == 1 else 0.0
                self.b3[l, m] = np.sqrt((l - m + 1) * (l - m + 2) * (2 * l + 2) * (1 + delta_1m) / (2 * l + 3))
                # b4
                self.b4[l, m] = np.sqrt((l - m + 1) * (l + m + 1) * (2 * l + 1) / (2 * l + 3))

    def _recursive_EF_vectorized(self, r_vecs):
        """
        向量化计算 E 和 F 矩阵
        输入: r_vecs [N, 3]  (N 个时间点)
        输出: E, F [L+3, L+3, N]
        """
        N = r_vecs.shape[0]
        x = r_vecs[:, 0]
        y = r_vecs[:, 1]
        z = r_vecs[:, 2]

        r_sq = x ** 2 + y ** 2 + z ** 2
        r = np.sqrt(r_sq)

        dim = self.L + 3
        # E, F 增加一个维度 N
        E = np.zeros((dim, dim, N))
        F = np.zeros((dim, dim, N))

        # --- 初始值 (公式 14) ---
        # 利用广播机制，scalar / array -> array
        E[0, 0, :] = self.R_ref / r
        F[0, 0, :] = 0.0

        factor_11 = np.sqrt(3) * self.R_ref ** 2 / (r ** 3)
        E[1, 1, :] = factor_11 * x
        F[1, 1, :] = factor_11 * y

        # E[1, 0] = sqrt(3) * z * (Re/r^2) * E[0,0]
        # 注意: E[0,0] 已经是 Re/r
        E[1, 0, :] = np.sqrt(3) * (z * self.R_ref) / r ** 2 * E[0, 0, :]
        F[1, 0, :] = 0.0

        # --- 预计算常用的几何项，避免循环中重复计算 ---
        term_x_base = (x * self.R_ref / r_sq)
        term_y_base = (y * self.R_ref / r_sq)
        term_z_base = (z * self.R_ref / r_sq)
        term_r_base = (self.R_ref ** 2 / r_sq)

        # --- 递归循环 (l 无法向量化，必须循环，但内部对 N 向量化) ---
        for l in range(2, self.L + 2):
            # 1. Sectorial terms (l = m)
            m = l
            frac_sectorial = np.sqrt((2 * m + 1) / (2.0 * m))

            # 向量化计算
            E[l, m, :] = frac_sectorial * (term_x_base * E[m - 1, m - 1, :] - term_y_base * F[m - 1, m - 1, :])
            F[l, m, :] = frac_sectorial * (term_x_base * F[m - 1, m - 1, :] + term_y_base * E[m - 1, m - 1, :])

            # 2. Zonal and Tesseral terms (l != m)
            # 这里的 m 循环是常数级别的 (最大 40-96)，比起 N (17280) 很小
            # 可以通过矩阵操作进一步优化，但在 Python 中 explicit loop over m is acceptable if N is vectorized
            for m in range(l):
                factor1 = np.sqrt((2 * l + 1) * (2 * l - 1) / ((l - m) * (l + m)))
                factor2 = np.sqrt((2 * l + 1) * (l - m - 1) * (l + m - 1) / ((l - m) * (l + m) * (2 * l - 3)))

                term1_E = factor1 * term_z_base * E[l - 1, m, :]
                term1_F = factor1 * term_z_base * F[l - 1, m, :]

                term2_E = factor2 * term_r_base * E[l - 2, m, :]
                term2_F = factor2 * term_r_base * F[l - 2, m, :]

                E[l, m, :] = term1_E - term2_E
                F[l, m, :] = term1_F - term2_F

        return E, F

    def compute_acceleration_vectorized(self, sat_pos_array, dC, dS):
        """
        计算一批卫星位置的加速度
        sat_pos_array: [N, 3]
        dC, dS: [L+1, L+1]
        返回: [N, 3] 加速度向量
        """
        N = sat_pos_array.shape[0]
        E, F = self._recursive_EF_vectorized(sat_pos_array)  # E, F shape: [dim, dim, N]

        gx = np.zeros(N)
        gy = np.zeros(N)
        gz = np.zeros(N)

        # 常数因子
        GM_R2 = self.GM / (self.R_ref ** 2)
        GM_2R2 = self.GM / (2 * self.R_ref ** 2)

        # 循环 l 和 m 进行累加
        # 这里的计算是对 N 个点同时进行的
        for l in range(2, self.L + 1):
            for m in range(l + 1):
                clm = dC[l, m]
                slm = dS[l, m]

                # 如果系数极小，跳过以节省时间
                if abs(clm) < 1e-20 and abs(slm) < 1e-20:
                    continue

                # 取出对应的 b 系数 (标量)
                b1 = self.b1[l, m]
                b2 = self.b2[l, m]
                b3 = self.b3[l, m]
                b4 = self.b4[l, m]

                # 取出 E, F 切片，形状为 (N,)
                # 公式用到下标 l+1
                E_lp1_m = E[l + 1, m, :]
                F_lp1_m = F[l + 1, m, :]

                # Z 分量 (所有 m 通用)
                # dg_z = (GM/R^2) * b4 * (-E[l+1, m] * C - F[l+1, m] * S)
                dg_z = GM_R2 * b4 * (-E_lp1_m * clm - F_lp1_m * slm)
                gz += dg_z

                if m == 0:
                    # m=0 特殊情况
                    E_lp1_1 = E[l + 1, 1, :]
                    F_lp1_1 = F[l + 1, 1, :]

                    dg_x = GM_R2 * (-b1 * E_lp1_1 * clm)
                    dg_y = GM_R2 * (-b1 * F_lp1_1 * clm)
                else:
                    # m > 0
                    E_lp1_mp1 = E[l + 1, m + 1, :]
                    F_lp1_mp1 = F[l + 1, m + 1, :]
                    E_lp1_mm1 = E[l + 1, m - 1, :]
                    F_lp1_mm1 = F[l + 1, m - 1, :]

                    term_x = (b2 * (-E_lp1_mp1 * clm - F_lp1_mp1 * slm) +
                              b3 * (E_lp1_mm1 * clm + F_lp1_mm1 * slm))
                    dg_x = GM_2R2 * term_x

                    term_y = (b2 * (-F_lp1_mp1 * clm + E_lp1_mp1 * slm) -
                              b3 * (F_lp1_mm1 * clm - E_lp1_mm1 * slm))
                    dg_y = GM_2R2 * term_y

                gx += dg_x
                gy += dg_y

        return np.column_stack((gx, gy, gz))

    def calculate_lgd_batch(self, sat1_pos_array, sat2_pos_array, dC, dS):
        """
        批量计算 LGD
        输入: sat1, sat2 形状均为 [N, 3]
        """
        # 1. 批量计算加速度
        acc1 = self.compute_acceleration_vectorized(sat1_pos_array, dC, dS)
        acc2 = self.compute_acceleration_vectorized(sat2_pos_array, dC, dS)

        # 2. 批量计算视线方向向量 (LOS Vector)
        diff_pos = sat1_pos_array - sat2_pos_array
        # axis=1 求范数，保持维度 (N, 1) 以便广播
        dist = np.linalg.norm(diff_pos, axis=1, keepdims=True)
        e_los = diff_pos / dist

        # 3. 投影
        diff_acc = acc1 - acc2
        # np.sum(A * B, axis=1) 等同于逐行 dot product
        lgd = np.sum(diff_acc * e_los, axis=1)

        return lgd


# --- 文件读取辅助函数优化 ---

def read_shc_file_structure(filename, max_degree, skip_header_token, is_itsg=False):
    """通用的 SHC 文件读取器"""
    clm = np.zeros((max_degree + 1, max_degree + 1))
    slm = np.zeros((max_degree + 1, max_degree + 1))

    try:
        with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        start_idx = 0
        for i, line in enumerate(lines):
            if skip_header_token in line:
                start_idx = i + 1
                break

        for line in lines[start_idx:]:
            if not line.strip(): continue
            parts = line.split()

            # 简单的格式检查
            if is_itsg and line.startswith('gfc'):
                # ITSG format: key l m c s ...
                l, m = int(parts[1]), int(parts[2])
                c_idx, s_idx = 3, 4
            elif not is_itsg and (line.startswith('GRCOF2') or len(parts) >= 7):
                # GRACE format: key l m c s ...
                l, m = int(parts[1]), int(parts[2])
                c_idx, s_idx = 3, 4
            else:
                continue

            if l <= max_degree and m <= max_degree:
                clm[l, m] = float(parts[c_idx])
                slm[l, m] = float(parts[s_idx])

    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None, None

    return clm, slm


def read_grace_fo_shc_file(filename, max_degree=96):
    clm, slm = read_shc_file_structure(filename, max_degree, '# End of YAML header', is_itsg=False)
    return {'clm': clm, 'slm': slm}


def read_itsg_data(filename, max_degree=40):
    clm, slm = read_shc_file_structure(filename, max_degree, 'end_of_head', is_itsg=True)
    return {'clm': clm, 'slm': slm}


# ==========================================
# 主程序
# ==========================================
if __name__ == "__main__":
    print("--- 开始 LGD 计算 ---")

    # =========================================================================
    # 1. 统一配置区域 (请在此处修改所有设置)
    # =========================================================================
    # [路径配置]
    BASE_DIR = r'G:\GROOPS\PNAS2020Workspace'  # 工作空间根目录
    OUTPUT_DIR = os.path.join(os.getcwd(), 'results')  # 结果保存目录

    # 背景场文件 (相对于 BASE_DIR 或绝对路径)
    STATIC_FILE_REL = r'grace_products\GSM-2_2020122-2020152_GRFO_UTCSR_BB01_0603'

    # [时间范围]
    START_DATE = '2020-06-01'
    END_DATE = '2020-08-30'

    # [计算参数]
    MAX_DEGREE = 40  # 球谐系数阶数
    SAMPLING_INTERVAL = 5  # 轨道采样间隔 (秒)
    GM_EARTH = 3.986004418e14  # 引力常数
    R_EARTH = 6378137.0  # 参考半径

    # [区域筛选]
    FILTER_PASS = False  # 是否只计算特定区域的轨道 (True/False)
    REGION_CFG = {
        'lon': (88, 92),  # 经度范围
        'lat': (22, 26),  # 纬度范围
        'lat_limit': (-80.0, 80.0),  # 轨道纬度限制
        'dir': 'asc'  # 轨道方向 'asc' (升轨) 或 'desc' (降轨)
    }
    # =========================================================================

    # 2. 准备工作
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    static_file_path = os.path.join(BASE_DIR, STATIC_FILE_REL)

    # 生成日期列表
    start_dt = datetime.strptime(START_DATE, '%Y-%m-%d')
    end_dt = datetime.strptime(END_DATE, '%Y-%m-%d')
    date_list = [(start_dt + timedelta(days=x)).strftime('%Y-%m-%d')
                 for x in range((end_dt - start_dt).days + 1)]

    # 3. 初始化 (加载背景场 & 初始化计算器)
    print(f"初始化计算器 (L={MAX_DEGREE})...")
    calculator = LGDCalculator(GM_EARTH, R_EARTH, MAX_DEGREE)

    print(f"读取背景场: {os.path.basename(static_file_path)}")
    itsg_ref = read_grace_fo_shc_file(static_file_path, max_degree=MAX_DEGREE)

    if itsg_ref['clm'] is None:
        print("错误: 无法读取背景场文件，请检查路径。")
        exit()

    # 4. 按天循环处理
    for date_str in date_list:
        print(f"\nProcessing: {date_str}")

        # --- (A) 读取 ITSG 文件 ---
        itsg_file = os.path.join(BASE_DIR, 'itsg_dataset', date_str[:7],
                                 f'ITSG-Grace_operational_Kalman_n40_{date_str}.gfc')

        if not os.path.exists(itsg_file):
            print(f"  [跳过] ITSG文件不存在: {itsg_file}")
            continue

        itsg_data = read_itsg_data(itsg_file, max_degree=MAX_DEGREE)
        if itsg_data['clm'] is None: continue

        dC = itsg_data['clm'] - itsg_ref['clm']
        dS = itsg_data['slm'] - itsg_ref['slm']

        # --- (B) 读取轨道数据 ---
        try:
            # 加载笛卡尔坐标 (用于计算 LGD)
            loader = OrbitLoader(date_str=date_str, groops_workspace_dir=BASE_DIR)
            orb_c = loader.load_orbit_data('gnv1b', 'C', 'cartesian')
            orb_d = loader.load_orbit_data('gnv1b', 'D', 'cartesian')

            # 提取位置并按间隔采样
            pos_c = np.array([o.position for o in orb_c])[::SAMPLING_INTERVAL]
            pos_d = np.array([o.position for o in orb_d])[::SAMPLING_INTERVAL]
            timestamps = np.array([o.timestamp for o in orb_c])[::SAMPLING_INTERVAL]

            # --- (C) 区域筛选 (可选) ---
            if FILTER_PASS:
                loader_geo = OrbitLoader(date_str=date_str, groops_workspace_dir=BASE_DIR)
                orb_geo = loader_geo.load_orbit_data('gnv1b', 'C', 'geodetic')
                pos_geo = np.array([o.position for o in orb_geo])[::SAMPLING_INTERVAL]

                lonlat = pos_geo[:, 0:2]
                _, indices = filter_complete_tracks_passing_region(
                    lonlat, REGION_CFG['lon'], REGION_CFG['lat'],
                    lat_limit=REGION_CFG['lat_limit'], separate=False, direction=REGION_CFG['dir']
                )

                if len(indices) == 0:
                    print("  [提示] 筛选后无有效轨道段，跳过。")
                    continue

                # 应用筛选
                pos_c = pos_c[indices]
                pos_d = pos_d[indices]
                timestamps = timestamps[indices]

            # --- (D) 批量计算 LGD ---
            print(f"  计算点数: {len(timestamps)} ...")
            # 注意: 这里调用的是之前优化好的 batch 方法
            # 如果你用的是原始未优化的类，请改回 calculate_lgd 并在内部加循环
            lgd_values = calculator.calculate_lgd_batch(pos_c, pos_d, dC, dS) * 1e9     # nm/s^2

            # --- (E) 保存结果 ---
            time_nums = np.array([datetime_to_datenum(t) for t in timestamps])
            save_data = np.column_stack((time_nums, lgd_values)).astype(object)

            save_path = os.path.join(OUTPUT_DIR, f'time-lgd-{date_str}.mat')
            savemat(save_path, {'time_lgd': save_data})
            print(f"  已保存至: {os.path.basename(save_path)}")

        except Exception as e:
            print(f"  [Error] 处理日期 {date_str} 时出错: {e}")