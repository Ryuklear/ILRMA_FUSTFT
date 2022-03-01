clear;
close all;

% Set parameters
Param.supervised = false;  %true�ŋ��t�����@/false�Ńu���C���h��@
Total_times = 5; %���������̌J��Ԃ���(����)
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
% L = 8192; % window length in STFT [points]
% eta = L/8; % shift length in STFT [points]

L = 4096; % window length in STFT [points]
eta = L/2; % shift length in STFT [points]

%L = 8192;       %FUSTFT�̃t���[���V�t�g��
%eta = L/4;      %FUSTFT�̃t���[���V�t�g��

number_of_padded_zeros = L;
FUSTFTpara.real_valued_flag = 1;        %0���ƕ��f���l�A�P���Ǝ����l
FUSTFTpara.phase_flag = 0;              %0 = simple STFT, 1 = phase-aware STFT
FUSTFTpara.frequency_type = 1;
FUSTFTpara.fftSize = L;
FUSTFTpara.shiftSize = eta;

FUSTFTpara.type = 2; %1 or 2 (1:�k�������STFT, 2:�k���搶��STFT)
%%%% type 1 %%%%%
FUSTFTpara.w_name = 'hamming';
%%%% type 2 %%%%%
FUSTFTpara.frag = 1; %�k���搶��STFT�Ŏg�p
window_num = 1; %1:�n�j���O�� 2:�n�~���O�� 3:�u���b�N�}���� 4:�ł��؂�K�E�X�� 5:�T�C����
if window_num  == 1 %���쐬
    FUSTFTpara.analysis_window = (0.5 - 0.5*cos(2*pi*((0:(L-1))'+ 0.5)/L))/sqrt(L);%�n�j���O��
    FUSTFTpara.dual_window = create_dual_window(FUSTFTpara.analysis_window,eta);
elseif window_num  == 2
    FUSTFTpara.analysis_window = (0.54 - 0.46*cos(2*pi*((0:(L-1))'+ 0.5)/L))/sqrt(L);%�n�~���O��
    FUSTFTpara.dual_window = create_dual_window(FUSTFTpara.analysis_window,eta);
elseif window_num  == 3
    FUSTFTpara.analysis_window = (0.42 - 0.5*cos(2*pi*((0:(L-1))'+ 0.5)/L) +  0.08*cos(4*pi*((0:(L-1))'+0.5)/L))/sqrt(L);%�u���b�N�}����
    FUSTFTpara.dual_window = create_dual_window(FUSTFTpara.analysis_window,eta);
elseif window_num  == 4
     FUSTFTpara.analysis_window =  exp(-4*pi*((0:(L-1))' - (L-1)/2).^2/(L^2))/sqrt(L); %�ł��؂�K�E�X��
    %FUSTFTpara.analysis_window =  exp(-4*pi*((0:(L-1))' - (Lh-1)/2).^2/(L^2))/sqrt(L);%�ł��؂�K�E�X��
    FUSTFTpara.dual_window = create_dual_window(FUSTFTpara.analysis_window,eta);
else
    FUSTFTpara.analysis_window = 0.5*sin(pi*((0:(L-1))'+ 0.5)/L)/sqrt(L);%�T�C����
    FUSTFTpara.dual_window = create_dual_window(FUSTFTpara.analysis_window,eta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawConv = true; % true or false (true: plot cost function values in each iteration and show convergence behavior, false: faster and do not plot cost function values)
normalize = true; % true or false (true: apply normalization in each iteration of ILRMA to improve numerical stability, but the monotonic decrease of the cost function may be lost. false: do not apply normalization)
%%%%%%%%%%%%%%%��Ė@�p�����[�^�ݒ�%%%%%%%%%%%%%%%
Sparse.parameter = 0.05;  %�X�p�[�X�x�p�����[�^��(���t����E�Ȃ�)
Sparse.T_im = 4096;  %4096�ȍ~�̐M����0
Sparse.nu = 4096;%8192;    %�d��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%�����̐ݒ�(Param.ns�Ɠ��������L�q)%%%%%
audio(1).name = './sound_source/Midi_sound/music1/Stem__Bass.wav';
audio(2).name = './sound_source/Midi_sound/music1/Stem__Piano.wav';
% audio(3).name = './sound_source/Midi_sound/music1/Stem__Drum.wav';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Param.supervised
   L_type = 1;   % 1 or 2 (1:rank-1 learning method, 2:rank-r learning method)
   ext = '.wav'; % extention of learning sound file
   %�����Ɠ������ɋL�q
   % For sound source 1
   Spv(1).file = './sound_source/BasisLearning_sound/music1/Basis_Bass2(G2-D4#)/'; 
   Spv(1).rootname = 'Basis_Bass_';  % file name (�ԍ��̑O�܂�)
   % For sound source 2
   Spv(2).file = './sound_source/BasisLearning_sound/music1/Basis_piano_diff/';    
   Spv(2).rootname = 'Basis_Piano_';     
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Impulse response(����1���珇�ԂɋL�q)���������ӁI
%30,150
%For sound source 1
impulse_name(1,:) = './Impulse_response/E2A(0.3s)/imp030/imp030.17'; %h11 
impulse_name(2,:) = './Impulse_response/E2A(0.3s)/imp030/imp030.29'; %h21
impulse_name(3,:) = './Impulse_response/E2A(0.3s)/imp030/imp030.23'; %h31

%For sound source 2
impulse_name(4,:) = './Impulse_response/E2A(0.3s)/imp150/imp150.17'; %h12 
impulse_name(5,:) = './Impulse_response/E2A(0.3s)/imp150/imp150.29'; %h22
impulse_name(6,:) = './Impulse_response/E2A(0.3s)/imp150/imp150.23'; %h32

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
mix = mix_sig'; %���������F(samples x channel)

if abs(max(max(mix))) > 1 % check clipping
    error('Cliping detected while mixing.\n');
end

Eval_ILRMA = zeros(Param.ns, 3,Total_times); %�S�Ă̐��l���ʂ��L�^�����
Eval_ILRMA_SP = zeros(Param.ns, 3, Total_times);

%%%%%�r�[���t�H�[�~���O�O����%%%%%
theta_1 = 40;
theta_2 = -40;
v = 340;
d = 0.1698;

%FUSTFT
y_1 = FUSTFT(mix_sig(1,:)',FUSTFTpara.analysis_window,eta,FUSTFTpara.frequency_type,FUSTFTpara.real_valued_flag,FUSTFTpara.phase_flag);
y_2 = FUSTFT(mix_sig(2,:)',FUSTFTpara.analysis_window,eta,FUSTFTpara.frequency_type,FUSTFTpara.real_valued_flag,FUSTFTpara.phase_flag);

%STFT
%y_1 = STFT_basic(mix_sig(1,:)',STFT.analysis_window, eta, STFT.frag);
%y_2 = STFT_basic(mix_sig(2,:)',STFT.analysis_window, eta, STFT.frag);
%y_3 = STFT_basic(mix_sig(3,:)',STFT.analysis_window, eta, STFT.frag);
X(:,:,1) = y_1;
X(:,:,2) = y_2;

I = size(y_1,1);
J = size(y_1,2);
lam = zeros(I,1);
for i = 1:I
    if i == 1 %���g��1�̎��A�l�Ԃ̎��ɂ͉e���Ȃ��ƍl�����邽�߁Ai=2�Ɠ����ɂ���
       fi = i/L * Param.fsResample;
       lam(i,1) = v / fi;
    else
       fi = i-1/L * Param.fsResample;
       lam(i,1) = v / fi;
    end
end

%�X�e�A�����O�x�N�g���쐬
A = zeros(Param.mic,Param.ns,I);
W_opt2 = zeros(Param.ns,Param.mic,I);
for i = 1:I
    A(:,:,i) = eye(2);
% A(:,:,i) = [1,1; 
%     exp(-1i*(2*pi*d*sin(theta_1))/lam(i,1)),exp(-1i*(2*pi*d*sin(theta_2))/lam(i,1));
%     exp(-1i*(2*2*pi*d*sin(theta_1))/lam(i,1)),exp(-1i*(2*2*pi*d*sin(theta_2))/lam(i,1))];
%&W_opt2(:,:,i) = pinv(A(:,:,i));
W_opt2(:,:,i) = eye(2);
end

Y = zeros(Param.mic,1,J);
R = zeros(Param.mic,Param.mic,I);

for i = 1:I
    R_sum = zeros(Param.mic);
   for j = 1:J
      %Y(:,:,j) = [y_1(i,j)';y_2(i,j)';y_3(i,j)'];
      Y(:,:,j) = [y_1(i,j)';y_2(i,j)'];
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
%�\������
X = permute(X,[3,2,1]);
Y_re = zeros(I,J,Param.ns);
for i = 1:I
    for n = 1:Param.ns
       Y_re(i,:,n) = W_opt2(n,:,i) * X(:,:,i);
    end
end

Y_re_PB = projectionBack(Y_re,A_opt);
[FUSTFTpara.diagonal_components,FUSTFTpara.non_diagonal_components] = create_tridiagonal_systems(FUSTFTpara.analysis_window,eta,FUSTFTpara.frequency_type);
FUSTFTpara.periodicity_flag = 0;
Result_signal(:,1) = Inv_FUSTFT(Y_re_PB(:,:,1,1),FUSTFTpara.analysis_window,FUSTFTpara.diagonal_components,FUSTFTpara.non_diagonal_components,eta,Param.sig_length,FUSTFTpara.frequency_type,FUSTFTpara.real_valued_flag,FUSTFTpara.phase_flag,FUSTFTpara.periodicity_flag);
Result_signal(:,2) = Inv_FUSTFT(Y_re_PB(:,:,2,1),FUSTFTpara.analysis_window,FUSTFTpara.diagonal_components,FUSTFTpara.non_diagonal_components,eta,Param.sig_length,FUSTFTpara.frequency_type,FUSTFTpara.real_valued_flag,FUSTFTpara.phase_flag,FUSTFTpara.periodicity_flag);

%Result_signal(:,1) = ISTFT_basic(Y_re_PB(:,:,1,1),STFT.dual_window,eta,Param.sig_length,STFT.frag);
%Result_signal(:,2) = ISTFT_basic(Y_re_PB(:,:,2,1),STFT.dual_window,eta,Param.sig_length,STFT.frag);
audiowrite(sprintf('signal1.wav'),Result_signal(:,1),Param.fsResample);
audiowrite(sprintf('signal2.wav'),Result_signal(:,2),Param.fsResample);
%%%%%%%%%%%%%%%%%%%%

%���t����or�u���C���h��������
for i = 1:Total_times
   if Param.supervised
   % Supervised separation based on ILRMA and ILRMA+Sparse    
      %[estimate]=Supervised_separation(mix,Param,Sparse,STFT,Spv,ext,L_type,musicSS,drawConv,normalize);
      %STFT�̎��̋��t���艹������
      
      [estimate] = Supervised_separation_FUSTFT(mix,Param,Sparse,FUSTFTpara,Spv,ext,L_type,musicSS,drawConv,normalize);
      %FUSTFT�̎��̋��t���艹������
   else
   % Blind source separation based on ILRMA and ILRMA+Sparse
   %FUSTFT�̎��̃u���C���h��������
    [estimate] = BSS_FUSTFT(mix,Param,Sparse,FUSTFTpara,drawConv,normalize,musicSS);
    %[estimate] = BSS(mix,Param,Sparse,FUSTFTpara,drawConv,normalize,musicSS);
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

%%%%%%%%%%�����o��%%%%%%%%%%
%refMic�̍�������
audiowrite(sprintf('%s/ObservedMixture.wav', outputDir), mix(:,refMic), Param.fsResample); % observed signal
%�^�̊e�����Ɗe��@�ł̐��艹���o��(�Ō�ɕ����������ʂ��o��)
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
%%%%%%%%%%���l�]���o��%%%%%%%%%%
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
        analysis_window = (0.5 - 0.5*cos(2*pi*((0:(window_length-1))'+ 0.5)/window_length))/sqrt(window_length); % ���K���Ώ̃n�j���O��
        dual_window = create_dual_window(analysis_window,time_skip); % �o�Α����쐬
        
    elseif window_num == 2
        analysis_window = (0.54 - 0.46*cos(2*pi*((0:(window_length-1))'+ 0.5)/window_length))/sqrt(window_length); % ���K���Ώ̃n�~���O��
        dual_window = create_dual_window(analysis_window,time_skip); % �o�Α����쐬
        
    elseif window_num == 3
        analysis_window = (0.42 - 0.5*cos(2*pi*((0:(window_length-1))'+ 0.5)/window_length) +  0.08*cos(4*pi*((0:(window_length-1))'+0.5)/window_length))/sqrt(window_length); % ���K���Ώ̃u���b�N�}����
        dual_window = create_dual_window(analysis_window,time_skip); % �o�Α����쐬
        
    elseif window_num == 4
        analysis_window = exp(-4*pi*((0:(window_length-1))' - (window_length-1)/2).^2/(window_length^2))/sqrt(window_length); % ���K���ł��؂�Ώ̃K�E�X��
        dual_window = create_dual_window(analysis_window,time_skip); % �o�Α����쐬
        
    elseif window_num == 5
        analysis_window = sin(pi*((0:(window_length-1))' + 0.5)/window_length)/sqrt(window_length); % �T�C����
        dual_window = create_dual_window(analysis_window,time_skip); % �o�Α����쐬
    end    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%