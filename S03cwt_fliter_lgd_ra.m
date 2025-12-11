clear;clc;

%% 单日计算
date_str = '2021-08-02';
freq = 0.2;     % 采样频率
process_cwt_fliter_signal_data(date_str, 'lgd', fs=freq, show_plots=true, f_high=12e-3);


%% 批量
clear;clc;
freq = 0.2;     % 采样频率

start_datenum = datenum(2020, 6, 1);
end_datenum = datenum(2020, 8, 30);

date_nums = start_datenum:end_datenum;

date_list = cellstr(datestr(date_nums, 'yyyy-mm-dd'));

for i=1: length(date_list)
    current_date = date_list{i};
    
    try
        process_cwt_fliter_signal_data(current_date, "lgd", fs=freq, show_plots=false, f_high=12e-3);
        fprintf('%s 数据处理完成\n', current_date);
    catch ME
        fprintf('错误: 处理 %s 数据时失败: %s\n', current_date, ME.message);
        fprintf('跳过该日期，继续处理下一天...\n');
    end
    
    fprintf('----------------------------------------\n');
end


