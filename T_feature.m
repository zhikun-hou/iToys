% 提取时域标量特征
% 滑动窗口在整个函数外面再套一层循环即可

function [Feature] = T_feature(ECG_train,Config)
    % ECG_train是一个0/1的脉冲信号
    arguments
        ECG_train { mustBeVector,mustBeNonempty }
        Config.FeatureName { mustBeTextScalar,mustBeNonempty }
        Config.SampleRate {mustBeNonempty,mustBeInteger,mustBeGreaterThan(Config.SampleRate,30)} = [] % 心电采样率
    end

    dt = 1000/Config.SampleRate;
    T_train = (0:1:numel(ECG_train)-1 ) * dt;

    % 心电图中各个spike的时刻
    R_idx = find(ECG_train==1);
    R_time = T_train(R_idx);

    % 计算R-R interval的列表
    RR_interval = diff(R_time);
    % 两个R-R interval间隔时间的差异
    RR_diff = diff(RR_interval);
    

    switch(Config.FeatureName)
        case "HR" % 心率的单位是每分钟
            % 需要去除掉左右可能的padding，以此来提高心率估计的准确性
            T_beats = sum(RR_interval);
            N_beats = numel(RR_interval);
            Feature = ( N_beats / T_beats) * 60*1000;
        case "HRs" % 自己起的名字，用来表示以秒为单位的心率
            T_beats = sum(RR_interval);
            N_beats = numel(RR_interval);
            Feature = ( N_beats / T_beats)*1000;
        case "RMSSD"
            Feature = rms(RR_diff);
        case "SDRR"
            Feature = std(RR_interval);
        case "pNN50"
            Feature = sum(RR_diff>50) / numel(RR_interval);
        otherwise

    end


end

