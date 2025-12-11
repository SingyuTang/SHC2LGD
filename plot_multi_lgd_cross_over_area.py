import os
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
from S02compute_grace_lgd import OrbitLoader
from S05plot_lgd_ra_cwt_filter import filter_complete_tracks_passing_region
from scipy.io import savemat, loadmat

from plot_single_lgd import load_orbit


def _add_month_annotations(ax, date_list):
    """添加月份标注"""
    # 将日期字符串转换为datetime对象
    date_objs = [datetime.strptime(date, '%Y-%m-%d') for date in date_list]

    # 分组日期到月份
    months = {}
    for date_obj in date_objs:
        month_key = date_obj.strftime('%B %Y')
        if month_key not in months:
            months[month_key] = []
        months[month_key].append(date_obj)

    # 计算每个月份在横轴上的位置
    month_positions = {}
    for month, dates_in_month in months.items():
        indices = [date_list.index(d.strftime('%Y-%m-%d')) for d in dates_in_month]
        avg_index = np.mean(indices)
        month_positions[month] = avg_index

    # 在图上添加月份标注
    y_lim = ax.get_ylim()
    y_pos = y_lim[1] + (y_lim[1] - y_lim[0]) * 0.12
    for month, pos in month_positions.items():
        x_pos = pos * 5  # 因为每个日期偏移5个单位
        ax.text(x_pos, y_pos, month, ha='center', va='bottom',
                fontweight='bold', fontsize=10, color='darkblue')


def _add_date_ticks(ax, date_list):
    """添加日期刻度标记"""
    y_lim = ax.get_ylim()

    for i, date_str in enumerate(date_list):
        x_pos = i * 5  # 计算横坐标位置

        # 添加淡色垂直线
        ax.axvline(x=x_pos, ymin=0, ymax=1, color='lightgray',
                   linewidth=0.8, alpha=0.7, zorder=0)

        # 提取日期和添加后缀
        day = int(date_str.split('-')[2])
        if 4 <= day <= 20 or 24 <= day <= 30:
            suffix = 'th'
        else:
            suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(day % 10, 'th')

        # 添加日期文本
        y_pos = y_lim[1] + (y_lim[1] - y_lim[0]) * 0.06
        ax.text(x_pos, y_pos, f'{day}{suffix}',
                ha='center', va='bottom', fontsize=8, color='black')

def load_lat_lgd_list(date_str, lon_range, lat_range, lat_limit, orbit_direction):
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
    orbitc, orbitd = load_orbit(date_str=date_str, groops_workspace=r'G:\GROOPS\PNAS2020Workspace', coord_type='geodetic')  # cartesian
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

def plot_multi_lgd_cross_over_area():
    date_list = [
        '2020-06-04', '2020-06-10', '2020-06-15', '2020-06-21', '2020-06-26',
        '2020-07-02', '2020-07-07', '2020-07-13', '2020-07-18', '2020-07-24', '2020-07-29',
        '2020-08-04', '2020-08-09', '2020-08-15', '2020-08-20'
    ]
    lon_range = (88, 92)  # 纬度范围（度）
    lat_range = (22, 26)  # 经度范围（度）
    lat_limit = (-80.0, 80.0)  # 轨道延申纬度范围
    orbit_direction = 'asc'  # 轨道方向

    input_dir = os.path.join(os.getcwd(), 'results')
    output_dir = os.path.join(os.getcwd(), 'results')

    # 加载LGD和轨道数据
    lat_list = []
    ori_lgd_list = []
    cwt_lgd_list = []
    for date_str in date_list:
        print(f"Processing {date_str}...")
        lats, ori_lgd, cwt_lgd = load_lat_lgd_list(date_str, lon_range, lat_range, lat_limit, orbit_direction)

        lat_list.append(lats)
        ori_lgd_list.append(ori_lgd)
        cwt_lgd_list.append(cwt_lgd)
    print("🎉 数据加载完成！")

    # 绘图
    fig = plt.figure(figsize=(8, 5))
    gs = fig.add_gridspec(1, 1)
    ax1 = fig.add_subplot(gs[0])    # 原始信号
    # ax2 = fig.add_subplot(gs[1])    # 小波信号

    # 添加指定纬度范围背景色
    ax1.axhspan(lat_range[0], lat_range[1], facecolor='gray', alpha=0.15)

    # 绘制每个日期的ori LGD曲线
    for i, (lats, signals) in enumerate(zip(lat_list, ori_lgd_list)):
        if len(lats) > 0 and len(signals) > 0:
            offset_signals = signals + i * 5  # 每日数据偏移5个单位
            ax1.scatter(offset_signals, lats, s=1, label=date_list[i])

    # 设置坐标轴标签
    data_label = 'ori LGD'
    ax1.set_xlabel(f'{data_label} (nm/s²)', fontsize=12)
    ax1.set_ylabel('Latitude (deg)', fontsize=12)

    # 添加月份和日期标记
    _add_month_annotations(ax1, date_list)
    _add_date_ticks(ax1, date_list)

    # 调整y轴范围，为标注留出空间
    y_lim = ax1.get_ylim()
    ax1.set_ylim(y_lim[0], y_lim[1] + (y_lim[1] - y_lim[0]) * 0.15)

    plt.tight_layout()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_filename = f'{date_list[0]}_{date_list[-1]}_ITSG_LGD_crossing_over_area_{timestamp}.png'
    save_path = os.path.join(output_dir, output_filename)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"💾 图形已保存: {save_path}")

    plt.show()

if __name__ == '__main__':
    plot_multi_lgd_cross_over_area()



