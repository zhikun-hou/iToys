% 提取频域特征
% 传入R_train，然后提取心率，计算功率谱、算特征

function [Feature] = F_feature(ECG_train,Config)
    % ECG_train是一个0/1的脉冲信号
    arguments
        ECG_train { mustBeVector,mustBeNonempty }
        Config.FeatureName { mustBeTextScalar,mustBeNonempty }
        Config.SampleRate {mustBeNonempty,mustBeInteger,mustBeGreaterThan(Config.SampleRate,30)} = [] % 心电采样率
        Config.Window { mustBePositive } % 滑动窗口的大小，单位是秒
    end

    dt = 1000/Config.SampleRate;

    N_window = round(Config.Window*1000 / dt);
    N_total  = numel(ECG_train);
    N_points = N_total-N_window+1;

    if(N_window<1)
        error("Window太小了");
    end

    % 心率计算 =============================================================
    HR = zeros(1,N_points);

    for i=1:N_points
        Segment = ECG_train(i:i+N_window-1);
        HR(i) = T_feature(Segment,"FeatureName","HRs","SampleRate",Config.SampleRate);
    end

%     plot(HR);

    % 功率谱提取 ===========================================================

    [PSD,PSD_f] = pwelch(HR,Config.SampleRate);
    LF = sum(PSD(PSD_f> 0.15 & PSD_f<=0.4));
    HF = sum(PSD(PSD_f>=0.04 & PSD_f<0.15));

    switch(Config.FeatureName)
        case "LF"
            Feature = LF;
        case "HF"
            Feature = HF;
        case "LF/HF"
            Feature = LF/HF;
        case "log(LF/HF)"
        case "ln(LF/HF)"
        case "log[LF/HF]"
        case "ln[LF/HF]"
            Feature = log(LF/HF);
        otherwise

    end


end

