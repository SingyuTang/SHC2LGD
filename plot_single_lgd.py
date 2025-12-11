import glob
import os
from datetime import datetime, timedelta
from S02compute_grace_lgd import OrbitLoader
from matplotlib import pyplot as plt
from scipy.special import lpmn, factorial
import numpy as np
import math
import re
from datetime import datetime, timedelta
from S02compute_grace_lgd import OrbitLoader
from S05plot_lgd_ra_cwt_filter import filter_complete_tracks_passing_region
from scipy.io import savemat, loadmat


def load_orbit(date_str="2021-07-27", groops_workspace=r'G:\GROOPS\WRR2021WorkspaceLRI', coord_type='cartesian'):
    """
    Load GNV1B orbit data from file
    """
    orbit_loader = OrbitLoader(date_str=date_str, groops_workspace_dir=groops_workspace)
    orbitc = orbit_loader.load_orbit_data('gnv1b', 'C', coord_type)
    orbitd = orbit_loader.load_orbit_data('gnv1b', 'D', coord_type)
    return orbitc, orbitd

def load_lat_lgd_list(date_str, lon_range, lat_range, lat_limit, orbit_direction):
    """
    Load LGD data and filter latitudes from file and filter tracks passing the region of interest
    """

    input_dir = os.path.join(os.getcwd(), 'results')

    ori_filename = os.path.join(input_dir, f'time-lgd-{date_str}.mat')  # 原始数据文件路径
    cwt_filename = os.path.join(input_dir, f'cwt_time-lgd-{date_str}.mat')  # 小波重构数据文件路径
    ori_var_name = 'time_lgd'
    cwt_var_name = 'cwt_lgd'

    if not os.path.exists(ori_filename):
        raise FileNotFoundError(f"原始数据文件不存在: {ori_filename}")
    if not os.path.exists(cwt_filename):
        raise FileNotFoundError(f"小波滤波数据文件不存在: {cwt_filename}")

    # 加载数据
    ori_data = loadmat(ori_filename)[ori_var_name].astype(np.float64)
    cwt_data = loadmat(cwt_filename)

    cwt_time = cwt_data['time'].squeeze()  # 当日时间，累积秒，如5、10、15、20、...
    cwt_signal = cwt_data[cwt_var_name].squeeze()

    ori_time = cwt_time
    ori_signal = ori_data[:, 1]

    # 确保信号长度一致
    min_len = min(len(ori_signal), len(cwt_signal))
    ori_signal = ori_signal[:min_len]
    cwt_signal = cwt_signal[:min_len]
    cwt_time = cwt_time[:min_len]

    # 提取轨道数据
    interval = 5  # 时间间隔（单位：秒）
    groops_workspace = os.path.dirname(os.getcwd())
    orbitc, orbitd = load_orbit(date_str=date_str, coord_type='geodetic', groops_workspace=groops_workspace)  # cartesian
    posc_geo = np.squeeze(np.array([obj.position for obj in orbitc])[::interval])
    posd_geo = np.squeeze(np.array([obj.position for obj in orbitd])[::interval])
    timestamps = np.array(list(map(lambda obj: obj.timestamp, orbitc))[::interval])

    # 过滤轨道数据
    lonlat = np.column_stack([posc_geo[:, 0], posc_geo[:, 1]])
    tracks, indices = filter_complete_tracks_passing_region(
        data=lonlat, lon_range=lon_range, lat_range=lat_range, lat_limit=lat_limit, separate=False,
        direction=orbit_direction)
    posc_geo = np.squeeze(posc_geo[indices])
    posd_geo = np.squeeze(posd_geo[indices])
    timestamps = timestamps[indices]

    return posc_geo[:, 1], ori_signal[indices], cwt_signal[indices]


def plot_cwt_lgd(
        date_str, lon_range, lat_range,
        orbit_direction='asc', lat_limit=(-80, 80)
):
    """
    Plot original LGD and cwt LGD comparison
    """
    input_dir = os.path.join(os.getcwd(), 'results')
    output_dir = os.path.join(os.getcwd(), 'results')
    print(f"Processing {date_str}...")
    lats, ori_lgd, cwt_lgd = load_lat_lgd_list(date_str, lon_range, lat_range, lat_limit, orbit_direction)

    # 绘图
    fig = plt.figure(figsize=(8, 5))
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))

    # 绘制第（1，1）个子图，origianl LGD
    ax1.scatter(ori_lgd,  lats, color='blue', s=1)
    ax1.set_xlabel(f'LGD (nm/s^2)')
    ax1.set_ylabel('latitude (°)')
    ax1.set_title(f'original LGD')

    # 绘制第（2，1）个子图，cwt LGD
    ax2.scatter(cwt_lgd,  lats, color='blue', s=1)
    ax2.set_xlabel(f'LGD (nm/s^2)')
    ax2.set_ylabel('latitude (°)')
    ax2.set_title(f'cwt LGD')

    plt.suptitle(f'滤波前后的LGD对比图-{date_str}', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'lgd_comparison_{date_str}.png'),
                dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')

    plt.show()

if __name__ == '__main__':
    date_str = '2020-07-07'
    lon_range = (88, 92)  # 经度范围（度）
    lat_range = (22, 26)  # 纬度范围（度）
    lat_limit = (-80.0, 80.0)  # 轨道延申纬度范围
    orbit_direction = 'asc'  # 轨道方向
    plot_cwt_lgd(date_str, lon_range, lat_range, lat_limit=lat_limit, orbit_direction=orbit_direction)

