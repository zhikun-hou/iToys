
ecg = EKG1;
fs = 200;
gr = 1;

%% Outputs
% qrs_amp_raw : amplitude of R waves amplitudes
% qrs_i_raw : index of R waves
% delay : number of samples which the signal is delayed due to the filtering
%% ================= Now Part of BioSigKit ==================== %%
if ~isvector(ecg)
  error('ecg must be a row or column vector');
end
ecg = ecg(:); % vectorize

%% ======================= Initialize =============================== %
delay = 0;
skip = 0;    % becomes one when a T wave is detected
m_selected_RR = 0;
mean_RR = 0;
ser_back = 0; 

%% ============ Noise cancelation(Filtering)( 5-15 Hz) =============== %%
 f1=5;                                                                      % cuttoff low frequency to get rid of baseline wander
 f2=15;                                                                     % cuttoff frequency to discard high frequency noise
 Wn=[f1 f2]*2/fs;                                                           % cutt off based on fs
 N = 3;                                                                     % order of 3 less processing
 [a,b] = butter(N,Wn);                                                      % bandpass filtering
 ecg_h = filtfilt(a,b,ecg);
 ecg_h = ecg_h/ max( abs(ecg_h));
%% ==================== derivative filter ========================== %%
% ------ H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2)) --------- %
if fs ~= 200
 int_c = (5-1)/(fs*1/40);
 b = interp1(1:5,[1 2 0 -2 -1].*(1/8)*fs,1:int_c:5);
else
 b = [1 2 0 -2 -1].*(1/8)*fs;   
end

 ecg_d = filtfilt(b,1,ecg_h);
 ecg_d = ecg_d/max(ecg_d);
%% ========== Squaring nonlinearly enhance the dominant peaks ========== %%
 ecg_s = ecg_d.^2;
%% ============  Moving average ================== %%
%-------Y(nt) = (1/N)[x(nT-(N - 1)T)+ x(nT - (N - 2)T)+...+x(nT)]---------%
ecg_m = conv(ecg_s ,ones(1 ,round(0.150*fs))/round(0.150*fs));
delay = delay + round(0.150*fs)/2;
%% ===================== Fiducial Marks ============================== %% 
% Note : a minimum distance of 40 samples is considered between each R wave
% since in physiological point of view no RR wave can occur in less than
% 200 msec distance
[pks,locs] = findpeaks(ecg_m,'MINPEAKDISTANCE',round(0.2*fs));
%% =================== Initialize Some Other Parameters =============== %%
LLp = length(pks);
% ---------------- Stores QRS wrt Sig and Filtered Sig ------------------%
qrs_c = zeros(1,LLp);           % amplitude of R        滑动平均 R波的幅度
qrs_i = zeros(1,LLp);           % index                 滑动平均 R波的位置
qrs_i_raw = zeros(1,LLp);       % amplitude of R        滤波信号 R波的位置
qrs_amp_raw= zeros(1,LLp);      % Index                 滤波信号 R波的幅度
% ------------------- Noise Buffers ---------------------------------%
nois_c = zeros(1,LLp);
nois_i = zeros(1,LLp);
% ------------------- Buffers for Signal and Noise ----------------- %


%% initialize the training phase (2 seconds of the signal) to determine the THR_SIG and THR_NOISE
% 这四个用于滑动平均信号
THR_SIG = max(ecg_m(1:2*fs))*1/3;  % 0.25 of the max amplitude  % THRESHOLD I1 
THR_NOISE = mean(ecg_m(1:2*fs))*1/2; % 0.5 of the mean signal is considered to be noise   % THRESHOLD I2  
SIG_LEV= THR_SIG; % SNKI
NOISE_LEV = THR_NOISE; % NPKI

%% Initialize bandpath filter threshold(2 seconds of the bandpass signal)
% 这四个用于滤波信号 
% 这两个是阈值的动态评估

THR_SIG1 = max(ecg_h(1:2*fs))*1/3; % THRESHOLD F1                            % 0.25 of the max amplitude 
THR_NOISE1 = mean(ecg_h(1:2*fs))*1/2; % THRESHOLD F2
% 这两个是Peak的动态评估
SIG_LEV1 = THR_SIG1;                                                        % Signal level in Bandpassed filter
NOISE_LEV1 = THR_NOISE1;                                                    % Noise level in Bandpassed filter
%% ============ Thresholding and desicion rule ============= %%
% Beat_C似乎是用来计数的，当前是第几个心跳
Beat_C = 0;                                                                 % Raw Beats
Beat_C1 = 0;                                                                % Filtered Beats
Noise_Count = 0;                                                            % Noise Counter
for i = 1 : LLp  % 遍历每个候选的peak
   %% ===== locate the corresponding peak in the filtered signal === %%
    % findpeaks是在ech_m上算的，滑动平均使用了卷积，所以长度会超出原始信号

   % 去除开始150ms和最后150ms的peak
   % 卷积在两边各伸出75ms，为了不包含这段空白的信息，就需要从150ms开始
    if locs(i)-round(0.150*fs)>= 1 && locs(i)<= length(ecg_h)
        % locs是在ech_m上算的，减去后面那段是为了变换，取出从真实的起点直到locs(i)时刻的滤波信号值
        % 然后计算这段时间内的最大值及其索引
        [y_i,x_i] = max(ecg_h(locs(i)-round(0.150*fs):locs(i)));
    else
        % 对于超过那两个范围的
        if i == 1 % 如果是在开头150ms
            % 计算从开始到peak的最大值
            [y_i,x_i] = max(ecg_h(1:locs(i)));
            % 这个search_back不是指回溯，是特指峰发生在卷积开头padding的情况
            ser_back = 1;
        elseif locs(i)>= length(ecg_h) 
            % 如果是为末尾150ms
            % 计算从peak直到结尾的最大值
            [y_i,x_i] = max(ecg_h(locs(i)-round(0.150*fs):end));
        end       
    end       
  %% ================= update the heart_rate ==================== %% 
    % 需要实时维护前8个RR-interval的均值
    if Beat_C >= 9        
        % 计算所谓的AVERAGE 1
        diffRR = diff(qrs_i(Beat_C-8:Beat_C));                                   % calculate RR interval
        mean_RR = mean(diffRR);                                            % calculate the mean of 8 previous R waves interval
        % 计算当前的RR-interval
        comp =qrs_i(Beat_C)-qrs_i(Beat_C-1);                                     % latest RR
        % 如果落在接受域外
        if comp <= 0.92*mean_RR || comp >= 1.16*mean_RR
     % ------ lower down thresholds to detect better in MVI -------- %
                % 说明是心率异常，可能要放宽阈值来检测是否存在Q波
                THR_SIG = 0.5*(THR_SIG); % I1
                THR_SIG1 = 0.5*(THR_SIG1); % F1      
        else
            % 落在接受域内，正常窦性心率
            m_selected_RR = mean_RR;                                     % The latest regular beats mean
        end 
          
    end
    
 %% == calculate the mean last 8 R waves to ensure that QRS is not ==== %%
    % test_m即所谓的AVERAGE 2
    % ？？？？？？？？？？？？根本没有按论文的方式计算AVERAGE 2
   if m_selected_RR % 是正常窦性心率
       test_m = m_selected_RR;                                         %if the regular RR availabe use it   
   elseif mean_RR && m_selected_RR == 0 % 不正常的情况，使用AVERAGE 1作为心率        为什么？？？？？？
       test_m = mean_RR;   
   else
       test_m = 0; % 似乎是心跳监护仪断开连接的情况，导致AVERAGE 1也为0
   end
        
    if test_m % 如果设备没有断开
        % 超过了拒绝域，说明有R波被忽略了，我们要在现在的R波和上一个R波中间插一个R波进去

        % 插一个R波进去，仍然有可能落在接受域外啊....那样的话就不应该算到AVERAGE 2里
          if (locs(i) - qrs_i(Beat_C)) >= round(1.66*test_m)                  % it shows a QRS is missed
              % 回溯这段时间的最大值，将其作为疑似的R波
              [pks_temp,locs_temp] = max(ecg_m(qrs_i(Beat_C)+ round(0.200*fs):locs(i)-round(0.200*fs))); % search back and locate the max in this interval
              locs_temp = qrs_i(Beat_C)+ round(0.200*fs) + locs_temp -1;      % location 
             
              if pks_temp > THR_NOISE % I2
                  % 超过了I2，说明确实存在一个漏掉的R波
                   Beat_C = Beat_C + 1;
                   qrs_c(Beat_C) = pks_temp;
                   qrs_i(Beat_C) = locs_temp;  
                   % 对阈值的调整在后面的部分
                  % ------------- Locate in Filtered Sig ------------- %
                   if locs_temp <= length(ecg_h)
                      [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.150*fs):locs_temp));
                   else
                      [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.150*fs):end));
                   end
                  % ----------- Band pass Sig Threshold ------------------%
                   if y_i_t > THR_NOISE1  % F2
                      Beat_C1 = Beat_C1 + 1;
                      qrs_i_raw(Beat_C1) = locs_temp-round(0.150*fs)+ (x_i_t - 1);% save index of bandpass 
                      qrs_amp_raw(Beat_C1) = y_i_t;                               % save amplitude of bandpass 
                      SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1;                      % when found with the second thres 
                   end
                   % not_noise = 1
                   % 论文里写的是用overall peak来更新，不应该是全体峰值吗？
                   SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV ;                       % when found with the second threshold             
             end             
          else
              % not_noise = 0
              % 说明中间没有遗漏心跳
          end
    end
  
    % 为什么test_m的end不一直延伸到最后？？？？？？


    %% ===================  find noise and QRS peaks ================== %%
    % 超过I1
    % 否则，如果小于I1超过I2  ->  更新对noise level的I/F评估
    %                           但是别的什么也没做？？？？可能逻辑在其他地方
    % 否则，如果小于I2        -> 多了一个noise峰，记录之，并更新对noise level的I/F评估

    if pks(i) >= THR_SIG % I1
      % ------ if No QRS in 360ms of the previous QRS See if T wave ------%
       if Beat_C >= 3
           % 间隔在360ms内，需要判断是不是T波误检为R波
          if (locs(i)-qrs_i(Beat_C)) <= round(0.3600*fs)
              Slope1 = mean(diff(ecg_m(locs(i)-round(0.075*fs):locs(i))));       % mean slope of the waveform at that position
              Slope2 = mean(diff(ecg_m(qrs_i(Beat_C)-round(0.075*fs):qrs_i(Beat_C)))); % mean slope of previous R wave
              if abs(Slope1) <= abs(0.5*(Slope2))                              % slope less then 0.5 of previous R
                 Noise_Count = Noise_Count + 1;
                 nois_c(Noise_Count) = pks(i);
                 nois_i(Noise_Count) = locs(i);
                 skip = 1;                                                 % T wave identification
                 % ----- adjust noise levels ------ %
                 NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
                 NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV; 
              else
                 skip = 0;
              end
            
           end
        end
        %---------- skip is 1 when a T wave is detected -------------- %
        if skip == 0    
          Beat_C = 1 + 1; % ？？？？？
          qrs_c(Beat_C) = pks(i);
          qrs_i(Beat_C) = locs(i);
        
        %--------------- bandpass filter check threshold --------------- %
          if y_i >= THR_SIG1  % F1
              Beat_C1 = Beat_C1 + 1; 
              if ser_back % 这个search_back不是指回溯，是特指峰发生在卷积开头padding的情况
                 qrs_i_raw(Beat_C1) = x_i;                                 % save index of bandpass 
              else
                 qrs_i_raw(Beat_C1)= locs(i)-round(0.150*fs)+ (x_i - 1);   % save index of bandpass 
              end
              qrs_amp_raw(Beat_C1) =  y_i;                                 % save amplitude of bandpass 
              SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;                       % adjust threshold for bandpass filtered sig
          end
         SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV ;                          % adjust Signal level
        end
              
    elseif (THR_NOISE <= pks(i)) && (pks(i) < THR_SIG)  % I2 I1
         % 说明是在超时中被诊断出来的漏诊峰
         NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;                        % adjust Noise level in filtered sig
         NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;                       % adjust Noise level in MVI       
    elseif pks(i) < THR_NOISE % I2
        Noise_Count = Noise_Count + 1; % 确实是一个noise峰
        nois_c(Noise_Count) = pks(i); % 记录这个noise峰在原始事件序列中的位置
        nois_i(Noise_Count) = locs(i); % 记录这个noise峰的幅度
        % 更新对于noise level的I/F估计
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;                         % noise level in filtered signal    
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;                        % adjust Noise level in MVI     
    end
               
    %% ================== adjust the threshold with SNR ============= %%
    % 正常情况下的浮动

    % 为什么要排除=0的时候？？？？？
    if NOISE_LEV ~= 0 || SIG_LEV ~= 0
        THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV)); % 更新I1
        THR_NOISE = 0.5*(THR_SIG); % 更新I2
    end
    
    %------ adjust the threshold with SNR for bandpassed signal -------- %
    if NOISE_LEV1 ~= 0 || SIG_LEV1 ~= 0
        THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));% 更新F1
        THR_NOISE1 = 0.5*(THR_SIG1); % 更新F2，追随I1浮动更新
    end
    
    
% ----------------------- reset parameters -------------------------- % 
    

    skip = 0;       
    ser_back = 0;

end
%% ======================= Adjust Lengths ============================ %%
qrs_i_raw = qrs_i_raw(1:Beat_C1);
qrs_amp_raw = qrs_amp_raw(1:Beat_C1);

qrs_c = qrs_c(1:Beat_C);
qrs_i = qrs_i(1:Beat_C);








