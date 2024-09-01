% 基于原始心电信号，输出脉冲图、滑动窗口心率图、RR间隔

% 代码修改子Pan Tompkin的原始论文
% 有很多Magic Number，有些似乎是生理约束，有些作者自己也没讲
% 【原始论文】Pan, Jiapu; Tompkins, Willis J. (March 1985). "A Real-Time QRS Detection Algorithm". IEEE Transactions on Biomedical Engineering. BME-32 (3): 230–236. doi:10.1109/TBME.1985.325532. PMID 3997178
% 【参考工具包】Complete Pan Tomkins Implementation ECG QRS detector  <-代码写得一坨

% 工具包写得一坨，比如：
% 没有按照论文中的方式实现AVERAGE 2
% 系数是有问题的，自己注释里是0.25，代码里就变成了1/3
% 它直接在findpeaks里使用了200ms的间隔，而不是按论文所说，对"确实无误的R波"之后使用200ms间隔
% 原始论文没说learning phase的具体实现，工具包直接用前两秒算了个值，根本没有learning
    % 而且工具包是直接findpeaks然后遍历找到的peak，原始论文的逻辑应该是遍历各个时刻的数据点
    % 因为原始论文是实时算法
% （29）上面的部分说，上一个R波到1.66*AVERAGE2之间的最大峰当作candidate，不是最大值，也不是直接选中

% 原始论文也充满歧义，比如：
% 在T波误检判别时，判别斜率是否为前一个期间波形的1/2，但没说是原始、滤波、滑动平均信号的哪一个

function [R_train,R_idx,R_amp] = PT(ECG_raw,Config)
    arguments
        ECG_raw { mustBeVector,mustBeNonempty  } % 原始心电数据，一维向量
        Config.Order {mustBeInteger,mustBeNonnegative} = 3 % FIR带通滤波器的阶数
        Config.SampleRate {mustBeNonempty,mustBeInteger,mustBeGreaterThan(Config.SampleRate,30)} = [] % 心电采样率
    end
    
    F_sample = Config.SampleRate;
    F_order  = Config.Order;

    % 带通滤波 ============================================================
    % 原始论文中，由于是实时监测，使用了因果律波器，一个低通一个高通级联
    % 这里既然是离线分析，就没必要了，直接FIR即可
    
    F_low = 5;
    F_high = 15; % 约束了SampleRate不能小于等于30
    
    F_normalize = [F_low,F_high]*2 / F_sample;
    [F_a,F_b] = butter(F_order,F_normalize);

    ECG_filtered = filtfilt(F_a,F_b,ECG_raw);

    % 差分信号 =============================================================

    D_b = [1 2 0 -2 -1].*(1/8)*F_sample;  
    ECG_derivate = filtfilt(D_b,1,ECG_filtered);
    
    % 整合信号 =============================================================
    
    % 工具包里是直接使用卷积实现的滑动平均，我直接用了movmean函数
    ECG_square = ECG_derivate.^2;
    ECG_avg = movmean(ECG_square,round(0.15*F_sample));
    
    % 参数部分 =============================================================
    
    Args = {};

    % 心电的两个样本点之间有多少ms的间隔
    dt = 1000/F_sample;

    % 两个R波之间小于200ms -> 不可能，将当前锋判定为噪声
    Args.P2P_noise = ceil(200/dt);

    % 两个R波之间小于360ms -> 有可能是R波也可能是T波，进行【误诊判断】
    Args.P2P_T = floor(360/dt);

    % 超过1.66倍AVERAGE 2时长无峰 -> 可能存在【漏检】，将上一个R到当前时刻间的最大峰作为
    Args.AVERAGE_1 = []; % 当前心率的平均R-R interval，间接衡量即时心率
    Args.AVERAGE_2 = []; % 正常心率的平均R-R interval，间接衡量平稳心率
    Args.RR_low  = 0.92; % 异常心率拒绝域的边界
    Args.RR_high = 1.16; % 同上
    Args.RR_miss = 1.66; % 心跳漏检的边界

    % 用于幅度判别的参数 ======
    % F表示滤波后的信号，I表示求导、平方、滑动平均之后的信号
    Args.SPKI = 0;
    Args.SPKF = 0;
    Args.NPKI = 0;
    Args.NPKF = 0;
    Args.THRESHOLD_I1 = 0;
    Args.THRESHOLD_I2 = 0;
    Args.THRESHOLD_F1 = 0;
    Args.THRESHOLD_F2 = 0;
    % 学习阶段
    Args.N_learn = round(2000 / dt); % 2s

    % 数据部分 =============================================================

    Context = {};

    % 所有R波的信息 ============
    Context.R_idx = zeros(1,0); % 所有R波的索引（滤波后、平方前）
    Context.R_amp = zeros(1,0); % 所有R波的幅度（滤波后、平方前）

    % 上一个R波的信息 ===========
    Context.Idx_last = [];  % 实验刚开始，还没有已经确定的R波索引(同上)，需要不断学习
    Context.Slope_last = []; % 存储上一个R波的最大斜率，用于T波误检的判别

    % 上一个R波至今所有峰的信息 =====
    Context.Cache_idx = zeros(1,0); % 由于可能漏检R，需要保留上一个R波至今的所有峰
    Context.Cache_amp = zeros(1,0);

    % 之前若干个RR-interval =======
    Context.RR_1 = zeros(1,0); % 前七个RR-interval，用于计算AVERAGE_1
    Context.RR_2 = zeros(1,0); % 前八个落在"正常域"内的RR-interval，用于计算AVERAGE_2

    % 开始检峰 =============================================================
    
    ECG_dir_last = 0; % 之前的方向，1代表上升，-1代表下降

    for i=2:numel(ECG_filtered) % 信号序列的原始长度，ECG_raw、ECG_filtered、ECG_derivate共用

        % R波漏诊鉴别 ======================================================

        % 将该模块放在最前
        % 1. 可以按照顺序，自然的将漏掉的R波插入到下一个R波之前
        % 2. 根据原论文的意思，一旦超过了166%的AVERAGE 2，就立刻寻找漏检的波形
        % 3. 如果不是按照2，而是在两个峰之前插入一个，有可能存在多个漏检

        % 确保AVERAGE_1和AVERAGE_2存在，拥有漏诊的判定标准
        if(not(isempty(Args.AVERAGE_2)))
            if(i-Context.Idx_last>=Args.RR_miss*Args.AVERAGE_2) % 超过时间阈值，可能存在漏检
                % 若漏诊，将上个峰到当前时刻之间最大的peak当作candidate
                [~,Cache_getter] = max(Context.Cache_amp);
                Context.Candidate_idx = Context.Cache_idx(Cache_getter);
                Context.CandidateF_amp = Context.Cache_amp(Cache_getter);
                Context.CandidateD_amp = ECG_derivate(Context.Candidate_idx);
                Context.CandidateI_amp = ECG_avg(Context.Candidate_idx);

                [Context,Args] = Step("past",Context,Args);
            end
        end
        
        % 极大值检测模块 ===================================================

        ECG_dir_new = sign(ECG_filtered(i)-ECG_filtered(i-1));

        if(ECG_dir_last==1 && ECG_dir_new<1) % 上升沿转平/下降沿，局部极大值出现
            is_local_maximum = true;
        else
            is_local_maximum = false;
        end
        % 更新当前信号的方向
        ECG_dir_last = ECG_dir_new;
        % 存在极值的话，需要提取出来并进行后续判别，不存在直接跳过该时刻
        if(is_local_maximum)
            % 注意：由于在i时刻检测到方向改变，i-1时刻才是真正的极大值点
            Context.Candidate_idx = i-1;
            Context.CandidateF_amp = ECG_filtered(Context.Candidate_idx);
            Context.CandidateD_amp = ECG_derivate(Context.Candidate_idx);
            Context.CandidateI_amp = ECG_avg(Context.Candidate_idx);

            [Context,Args] = Step("current",Context,Args);
        end
    
    end

    R_idx = Context.R_idx;
    R_amp = Context.R_amp;
    R_train = zeros(1,numel(ECG_raw));
    R_train(R_idx) = 1;

end

function [Context,Args] = Step(Candidate_from,Context,Args)

    % 基于时间判别 =====================================================
    if(Context.Candidate_idx<=Args.N_learn)
        % 最开始的学习阶段，不把任何峰判定为R波
        % 这段时间的峰只用来学习阈值
        Vote_time = -1;
    elseif(isempty(Context.Idx_last))
        % 学习阶段结束，并且之前一个Q波都没有
        Vote_time = 1;
    elseif(Context.Candidate_idx-Context.Idx_last<=Args.P2P_noise)
        % 当前候选峰与前一个被认定为R波的峰相差小于200ms，认定为噪声
        Vote_time = -1;
    else
        % 除了必须<=200ms外，没有额外的约束
        Vote_time = 1;
    end

    % 基于幅度判别 =====================================================
    
    % 注意：candidate可能来自于当前的峰，也有可能因为超时而来自过去的峰
    Above_F1 = false;
    Above_F2 = false;
    Above_I1 = false;
    Above_I2 = false;
    % 把比较的逻辑提出来，不能直接&&，也不能elseif
    % 因为PKI和PKF都要更新，所以必须全都比较完才行
    if(Context.CandidateF_amp>=Args.THRESHOLD_F1)
        Above_F1 = true;
    end
    if(Context.CandidateF_amp>=Args.THRESHOLD_F2)
        Above_F2 = true;
    end
    if(Context.CandidateI_amp>=Args.THRESHOLD_I1)
        Above_I1 = true;
    end
    if(Context.CandidateI_amp>=Args.THRESHOLD_I2)
        Above_I2 = true;
    end
    
    if(Candidate_from=="current")
        if(Above_F1 && Above_I1)
            % 两个都超过阈值1，毫无疑问是R波
            Vote_amp = 1;
        else
            % 没有超过阈值1，可能是R波可能是噪声，但姑且认为是噪声
            % 如果超时的话，会作为past自动返回来
            Vote_amp = -1;
        end
    elseif(Candidate_from=="past")
        if(Above_F2 && Above_I2)
            % 从幅度来看，之前似乎确实漏诊了一个
            % 但是，也有可能是T波，还需要经过一轮检验
            Vote_amp = 1;
        else
            % 之前没有漏诊，铁定不是R波
            Vote_amp = -1;
        end
    end
    
    % T波误诊鉴别 ======================================================
    
    if(isempty(Context.Idx_last))
        Vote_T = 1; % 第一个峰，不用考虑前一个峰的时间
    elseif(Context.Candidate_idx-Context.Idx_last<=Args.P2P_T)
        % 小于360ms，可能是T波误诊
        
        % 斜率比较算法，这里似乎应该使用ECG_derivate
        if(Context.CandidateD_amp<=0.5*Context.Slope_last)
            Vote_T = 1; % 斜率比较小，说明是T波而不是R波，误诊了
        else % 没有误诊，
            Vote_T = -1;
        end
    else
        % 大于360ms，不用担心是T波
        Vote_T = 1;
    end

    % 判别结果 =========================================================

    % 三轮判别，有任何一轮拒绝都不行
    is_R = Vote_time>0 && Vote_amp>0 && Vote_T>0;

    % 更新 ============================================================

    % 更新噪声/信号幅值的估计
    % 需要考虑一个问题：有可能F/I两个阈值，一个通过了一个没通过，而只有两个都通过才is_R
    % 所以要从is_R的花括号里提出来
    % 此外，学习阶段的峰我们不将其作为R波，只是将其用于更新阈值，所以也要与is_R解耦
    if(Candidate_from=="current")
        if(Above_F1)
            Args.SPKF = 0.125 * Context.CandidateF_amp + 0.875 * Args.SPKF;
        else
            Args.NPKF = 0.125 * Context.CandidateF_amp + 0.875 * Args.NPKF;
        end
        Args.THRESHOLD_F1 = Args.NPKF + 0.25 * (Args.SPKF - Args.NPKF);
        Args.THRESHOLD_F2 = 0.5 * Args.THRESHOLD_F1;
        if(Above_I1)
            Args.SPKI = 0.125 * Context.CandidateI_amp + 0.875 * Args.SPKI;
        else
            Args.NPKI = 0.125 * Context.CandidateI_amp + 0.875 * Args.NPKI;
        end
        Args.THRESHOLD_I1 = Args.NPKI + 0.25 * (Args.SPKI - Args.NPKI);
        Args.THRESHOLD_I2 = 0.5 * Args.THRESHOLD_I1;
    elseif(Candidate_from=="past") % 公式（21）
        if(Above_F2)
            Args.SPKF = 0.25 * Context.CandidateF_amp + 0.75 * Args.SPKF;
        else
            Args.NPKF = 0.25 * Context.CandidateF_amp + 0.75 * Args.NPKF;
        end
        if(Above_I2)
            Args.SPKI = 0.25 * Context.CandidateI_amp + 0.75 * Args.SPKI;
        else
            Args.NPKI = 0.25 * Context.CandidateI_amp + 0.75 * Args.NPKI;
        end
        if(Above_F2 && Above_I2) % 公式（22）（23）
            % 这里只需要更新THRESHOLD_1，因为算完past之后就算current，用的是THRESHOLD_1
            Args.THRESHOLD_I1 = 0.5 * Args.THRESHOLD_I1;
            Args.THRESHOLD_F1 = 0.5 * Args.THRESHOLD_F1;
        end
    end

    if(is_R)
        % 更新RR-interval列表
        % 根据论文的说法，使用前两个心跳构成的一个RR-interval就能进行初始化
        if(not(isempty(Context.Idx_last)))
            RR_current = Context.Candidate_idx-Context.Idx_last;
            
            if(numel(Context.RR_1)<8)
                % 还没有填满8个
                Context.RR_1(end+1) = RR_current;
                Context.RR_2(end+1) = RR_current;
    
                Args.AVERAGE_1 = mean(Context.RR_1);
                Args.AVERAGE_2 = mean(Context.RR_2);
            else
                % 环形缓冲区，即时心率
                Context.RR_1 = [Context.RR_1(2:end),RR_current];
                Args.AVERAGE_1 = mean(Context.RR_1);
                
                % 平稳心率
                if(RR_current>=Args.RR_low*Args.AVERAGE_2 && RR_current<=Args.RR_high*Args.AVERAGE_2)
                    Context.RR_2 = [Context.RR_2(2:end),RR_current];
                    Args.AVERAGE_2 = mean(Context.RR_2);
                end
            end
        end

        % 将当前R波作为最后一个R波，用于200ms、360ms、漏检判别
        Context.Idx_last = Context.Candidate_idx;
        Context.Slope_last = Context.CandidateD_amp;
        % 把当前R波放进所有R波的列表里
        Context.R_idx(end+1) = Context.Candidate_idx;
        Context.R_amp(end+1) = Context.CandidateF_amp;
        % 清空缓存的、用于检漏的峰
        Context.Cache_idx = zeros(1,0);
        Context.Cache_amp = zeros(1,0);
    else
        % 不是R波，需要暂存起来，以便漏检时取用
        Context.Cache_idx(end+1) = Context.Candidate_idx;
        Context.Cache_amp(end+1) = Context.CandidateF_amp;
    end
    

end


