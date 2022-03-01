clear;
close all;

% Set parameters
Param.supervised = false;  %trueで教師あり手法
Total_times = 10; %分離実験の繰り返し回数(平均)
%seed = 1; % pseudo random seed
refMic = 1; % reference microphone for estimated signal
Param.seconds = 13; % length of music
Param.fsResample = 16000; % resampling frequency [Hz]
Param.ns = 2; % number of sources
Param.mic = 2; %number of microphones
Param.nb = 30; % number of bases (for type=1, nb is # of bases for "each" source. for type=2, nb is # of bases for "all" sources)
Param.it = 100; % number of iterations (define by checking convergence behavior with drawConv=true)
Param.type = 1; % 1 or 2 (1: ILRMA w/o partitioning function, 2: ILRMA with partitioning function)
%%%%%%%%%%%%%%%STFT parameters%%%%%%%%%%%%%%%
L = 8192; % window length in STFT [points]
eta = L/4; % shift length in STFT [points]

%L = 8192;       %FUSTFTのフレームシフト数
%eta = L/4;      %FUSTFTのフレームシフト量


number_of_padded_zeros = L;
real_valued_flag = 1;        %0だと複素数値、１だと実数値
phase_flag = 0;              %0 = simple STFT, 1 = phase-aware STFT
analysis_window = (0.5 - 0.5*cos(2*pi*((0:(L-1))'+ 0.5)/L))/sqrt(L); % 正規化対称ハン窓



FUSTFTpara.type = 2; %1 or 2 (1:北村さんのSTFT, 2:北原先生のSTFT)
%%%% type 1 %%%%%
FUSTFTpara.w_name = 'hamming';
%%%% type 2 %%%%%
FUSTFTpara.frag = 1; %北原先生のSTFTで使用
window_num = 2; %1:ハニング窓 2:ハミング窓 3:ブラックマン窓 4:打ち切りガウス窓 5:サイン窓
if window_num  == 1 %窓作成
    FUSTFTpara.analysis_window = (0.5 - 0.5*cos(2*pi*((0:(L-1))'+ 0.5)/L))/sqrt(L);
    FUSTFTpara.dual_window = create_dual_window(FUSTFTpara.analysis_window,eta);
elseif window_num  == 2
    FUSTFTpara.analysis_window = (0.54 - 0.46*cos(2*pi*((0:(L-1))'+ 0.5)/L))/sqrt(L);
    FUSTFTpara.dual_window = create_dual_window(FUSTFTpara.analysis_window,eta);
elseif window_num  == 3
    FUSTFTpara.analysis_window = (0.42 - 0.5*cos(2*pi*((0:(L-1))'+ 0.5)/L) +  0.08*cos(4*pi*((0:(L-1))'+0.5)/L))/sqrt(L);
    FUSTFTpara.dual_window = create_dual_window(FUSTFTpara.analysis_window,eta);
elseif window_num  == 4
    FUSTFTpara.analysis_window = 1;
    FUSTFTpara.dual_window = create_dual_window(FUSTFTpara.analysis_window,eta);
else
    FUSTFTpara.analysis_window = 0.5*sin(pi*((0:(L-1))'+ 0.5)/L)/sqrt(L);
    FUSTFTpara.dual_window = create_dual_window(FUSTFTpara.analysis_window,eta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawConv = true; % true or false (true: plot cost function values in each iteration and show convergence behavior, false: faster and do not plot cost function values)
normalize = true; % true or false (true: apply normalization in each iteration of ILRMA to improve numerical stability, but the monotonic decrease of the cost function may be lost. false: do not apply normalization)
%%%%%%%%%%%%%%%提案法パラメータ設定%%%%%%%%%%%%%%%
Sparse.parameter = 0.05;  %スパース度パラメータλ(教師あり・なし)
Sparse.T_im = 4096;  %4096以降の信号が0
Sparse.nu = 8192;    %重み
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%音源の設定(Param.nsと同じだけ記述)%%%%%
audio(1).name = './sound_source/Midi_sound/music1/Cello.wav';
audio(2).name = './sound_source/Midi_sound/music1/Violin.wav';
% audio(3).name = './sound_source/Midi_sound/music1/Stem__Drum.wav';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Param.supervised
   L_type = 1;   % 1 or 2 (1:rank-1 learning method, 2:rank-r learning method)
   ext = '.wav'; % extention of learning sound file
   %音源と同じ順に記述
   % For sound source 1
   Spv(1).file = './sound_source/BasisLearning_sound/music1/Basis_Bass2(G2-D4#)/'; 
   Spv(1).rootname = 'Basis_Bass_';  % file name (番号の前まで)
   % For sound source 2
   Spv(2).file = './sound_source/BasisLearning_sound/music1/Basis_piano_diff/';    
   Spv(2).rootname = 'Basis_Piano_';     
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Impulse response(音源1から順番に記述)書き方注意！
%For sound source 1
impulse_name(1,:) = './Impulse_response/E2A(0.3s)/imp050/imp050.17'; %h11 
impulse_name(2,:) = './Impulse_response/E2A(0.3s)/imp050/imp050.29'; %h21
impulse_name(3,:) = './Impulse_response/E2A(0.3s)/imp050/imp050.23'; %h31

%For sound source 2
impulse_name(4,:) = './Impulse_response/E2A(0.3s)/imp130/imp130.17'; %h12 
impulse_name(5,:) = './Impulse_response/E2A(0.3s)/imp130/imp130.29'; %h22
impulse_name(6,:) = './Impulse_response/E2A(0.3s)/imp130/imp130.23'; %h32

%For sound source 3
% impulse_name(7,:) = './Impulse_response/E2A(0.3s)/imp090/imp090.17'; %h13
% impulse_name(8,:) = './Impulse_response/E2A(0.3s)/imp090/imp090.30'; %h23
% impulse_name(9,:) = './Impulse_response/E2A(0.3s)/imp090/imp090.23'; %h33
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix random seed

%RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed))

%Read inpulse responses
[audio] = reading_impulse(audio,impulse_name,Param,L);
%Make mixed sounds

[mix_sig,audio,Param] = audio_rc(audio,Param,FUSTFTpara);
%[mix_sig,audio,Param] = audio_rc(audio,Param,STFT);
musicSS = zeros(Param.ns,length(audio(1).resample)); %(ns x samples)

for n = 1:Param.ns
musicSS(n,:) = audio(n).resample';
end
mix = mix_sig'; %混合音源：(samples x channel)

if abs(max(max(mix))) > 1 % check clipping
    error('Cliping detected while mixing.\n');
end

Eval_ILRMA = zeros(Param.ns, 3,Total_times); %全ての数値結果が記録される
Eval_ILRMA_SP = zeros(Param.ns, 3, Total_times);

%%%%%ビームフォーミング前処理%%%%%
theta_1 = 40;
theta_2 = -40;
v = 340;
d = 0.1698;

%FUSTFT
frequency_type =1;
y_1 = FUSTFT(mix_sig(1,:)',analysis_window,eta,frequency_type,real_valued_flag,phase_flag)
y_2 = FUSTFT(mix_sig(1,:)',analysis_window,eta,frequency_type,real_valued_flag,phase_flag)
y_3 = FUSTFT(mix_sig(1,:)',analysis_window,eta,frequency_type,real_valued_flag,phase_flag)




%STFT
%y_1 = STFT_basic(mix_sig(1,:)',STFT.analysis_window, eta, STFT.frag);
%y_2 = STFT_basic(mix_sig(2,:)',STFT.analysis_window, eta, STFT.frag);
%y_3 = STFT_basic(mix_sig(3,:)',STFT.analysis_window, eta, STFT.frag);
X(:,:,1) = y_1;
X(:,:,2) = y_2;
X(:,:,3) = y_3;
I = size(y_1,1);
J = size(y_1,2);
lam = zeros(I,1);
for i = 1:I
    if i == 1 %周波数1の時、人間の耳には影響ないと考えられるため、i=2と同じにする
       fi = i/L * Param.fsResample;
       lam(i,1) = v / fi;
    else
       fi = i-1/L * Param.fsResample;
       lam(i,1) = v / fi;
    end
end

%ステアリングベクトル作成
A = zeros(Param.mic,Param.ns,I);
W_opt2 = zeros(Param.ns,Param.mic,I);
for i = 1:I
A(:,:,i) = [1,1; 
    exp(-1i*(2*pi*d*sin(theta_1))/lam(i,1)),exp(-1i*(2*pi*d*sin(theta_2))/lam(i,1));
    exp(-1i*(2*2*pi*d*sin(theta_1))/lam(i,1)),exp(-1i*(2*2*pi*d*sin(theta_2))/lam(i,1))];
W_opt2(:,:,i) = pinv(A(:,:,i));
end

Y = zeros(Param.mic,1,J);
R = zeros(Param.mic,Param.mic,I);

for i = 1:I
    R_sum = zeros(Param.mic);
   for j = 1:J
      Y(:,:,j) = [y_1(i,j)';y_2(i,j)';y_3(i,j)'];
      R_sum = R_sum + Y(:,:,j)*Y(:,:,j)';
   end
   R(:,:,i) = R_sum /J;
end

Param.W_opt = zeros(Param.ns,Param.mic,I);
A_opt = zeros(Param.mic,Param.ns,I);
for i = 1:I
W = pinv(R(:,:,i))*A(:,:,i)*pinv((A(:,:,i)'*pinv(R(:,:,i))*A(:,:,i)))*eye(Param.ns);
Param.W_opt(:,:,i) = W';
A_opt(:,:,i) = pinv(Param.W_opt(:,:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%予備実験
X = permute(X,[3,2,1]);
Y_re = zeros(I,J,Param.ns);
for i = 1:I
    for n = 1:Param.ns
       Y_re(i,:,n) = W_opt2(n,:,i) * X(:,:,i);
    end
end

Y_re_PB = projectionBack(Y_re,A_opt);
Result_signal(:,1) = ISTFT_basic(Y_re_PB(:,:,1,1),STFT.dual_window,eta,Param.sig_length,STFT.frag);
Result_signal(:,2) = ISTFT_basic(Y_re_PB(:,:,2,1),STFT.dual_window,eta,Param.sig_length,STFT.frag);
audiowrite(sprintf('signal1.wav'),Result_signal(:,1),Param.fsResample);
audiowrite(sprintf('signal2.wav'),Result_signal(:,2),Param.fsResample);
%%%%%%%%%%%%%%%%%%%%

for i = 1:Total_times
   if Param.supervised
   % Supervised separation based on ILRMA and ILRMA+Sparse    
      [estimate] = Supervised_separation(mix,Param,Sparse,STFT,Spv,ext,L_type,musicSS,drawConv,normalize);
   else
   % Blind source separation based on ILRMA and ILRMA+Sparse
      [estimate] = BSS(mix,Param,Sparse,STFT,drawConv,normalize,musicSS);
   end

%Evaluation
 ILRMA_temp = zeros(Param.ns,3,Param.mic);
 ILRMA_SP_temp = zeros(Param.ns,3,Param.mic);

   for m = 1:Param.mic
    %[SDR,SIR,SAR]
    [ILRMA_temp(:,1,m),ILRMA_temp(:,2,m),ILRMA_temp(:,3,m),~] = bss_eval_sources(estimate(m).T_signal_ILRMA',musicSS);
    [ILRMA_SP_temp(:,1,m),ILRMA_SP_temp(:,2,m),ILRMA_SP_temp(:,3,m),~] = bss_eval_sources(estimate(m).T_signal_ILRMA_SP',musicSS);
   end
   
  Eval_ILRMA(:,:,i) = sum(ILRMA_temp,3) / Param.mic; %[N x Evaluations x Total_times]
  Eval_ILRMA_SP(:,:,i) = sum(ILRMA_SP_temp,3) / Param.mic;
end

% Output separated signals
outputDir = sprintf('./output');
if ~isfolder( outputDir )
    mkdir( outputDir );
end

%%%%%%%%%%音源出力%%%%%%%%%%
%refMicの混合音源
audiowrite(sprintf('%s/ObservedMixture.wav', outputDir), mix(:,refMic), Param.fsResample); % observed signal
%真の各音源と各手法での推定音源出力(最後に分離した結果を出力)
for n = 1:Param.ns
Original_name = ['Original_Source',num2str(n),'.wav'];    
ILRMA_name = ['ES',num2str(n),'_ILRMA.wav'];
ILRMA_Sp_name = ['ES',num2str(n),'_ILRMA+SP.wav'];

audiowrite(sprintf('%s/%s', outputDir,Original_name), audio(n).resample, Param.fsResample);
audiowrite(sprintf('%s/%s', outputDir,ILRMA_name), estimate(refMic).T_signal_ILRMA(:,n), Param.fsResample); 
audiowrite(sprintf('%s/%s', outputDir,ILRMA_Sp_name), estimate(refMic).T_signal_ILRMA_SP(:,n), Param.fsResample); 
end
fprintf('The files are saved in "./output".\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%数値評価出力%%%%%%%%%%
Eval_ILRMA = permute(Eval_ILRMA,[1,3,2]); %[N x Total_times x Evaluations]
Eval_ILRMA_SP = permute(Eval_ILRMA_SP,[1,3,2]);
fprintf('ILRMA_SDR \n');
fprintf('%f\n',sum(Eval_ILRMA(:,:,1),2)/Total_times); fprintf('\n');
fprintf('ILRMA_SIR \n');
fprintf('%f\n',sum(Eval_ILRMA(:,:,2),2)/Total_times); fprintf('\n');
fprintf('ILRMA_SAR \n');
fprintf('%f\n',sum(Eval_ILRMA(:,:,3),2)/Total_times); fprintf('\n');
fprintf('ILRMA_SP_SDR \n');
fprintf('%f\n',sum(Eval_ILRMA_SP(:,:,1),2)/Total_times); fprintf('\n');
fprintf('ILRMA_SP_SIR \n');
fprintf('%f\n',sum(Eval_ILRMA_SP(:,:,2),2)/Total_times); fprintf('\n');
fprintf('ILRMA_SP_SAR \n');
fprintf('%f\n',sum(Eval_ILRMA_SP(:,:,3),2)/Total_times); fprintf('\n');

%%%%%%%%%%%%%%%Local function%%%%%%%%%%%%%%%
% Make window function in STFT_type2 
function[analysis_window,dual_window]  = Select_window(window_num,window_length,time_skip)

    if window_num == 1
        analysis_window = (0.5 - 0.5*cos(2*pi*((0:(window_length-1))'+ 0.5)/window_length))/sqrt(window_length); % 正規化対称ハニング窓
        dual_window = create_dual_window(analysis_window,time_skip); % 双対窓を作成
        
    elseif window_num == 2
        analysis_window = (0.54 - 0.46*cos(2*pi*((0:(window_length-1))'+ 0.5)/window_length))/sqrt(window_length); % 正規化対称ハミング窓
        dual_window = create_dual_window(analysis_window,time_skip); % 双対窓を作成
        
    elseif window_num == 3
        analysis_window = (0.42 - 0.5*cos(2*pi*((0:(window_length-1))'+ 0.5)/window_length) +  0.08*cos(4*pi*((0:(window_length-1))'+0.5)/window_length))/sqrt(window_length); % 正規化対称ブラックマン窓
        dual_window = create_dual_window(analysis_window,time_skip); % 双対窓を作成
        
    elseif window_num == 4
        analysis_window = exp(-4*pi*((0:(window_length-1))' - (window_length-1)/2).^2/(window_length^2))/sqrt(window_length); % 正規化打ち切り対称ガウス窓
        dual_window = create_dual_window(analysis_window,time_skip); % 双対窓を作成
        
    elseif window_num == 5
        analysis_window = sin(pi*((0:(window_length-1))' + 0.5)/window_length)/sqrt(window_length); % サイン窓
        dual_window = create_dual_window(analysis_window,time_skip); % 双対窓を作成
    end    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%