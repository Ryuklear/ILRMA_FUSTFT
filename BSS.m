function [estimate] = BSS(mix,Param,Sparse,STFT,drawConv,normalize,musicSS)
% [inputs]
%        mix: observed mixture (signal x channel)
%      Param: 各パラメータを含んだ構造体
%       type: without or with partitioning function (1: ILRMA without partitioning function (ILRMA1), 2: ILRMA with partitioning function (ILRMA2))
%     refMic: reference microphone for applying back projection
%   drawConv: plot cost function values in each iteration or not (true or false)
%  normalize: normalize variables in each iteration to avoid numerical divergence or not (true or false, normalization may collapse monotonic decrease of the cost function)
%
% [outputs]
%   estimate: 推定音源に関する構造体(estimate(M)まで存在)
%        >estimate(1).Y_ILRMA           : マイク1のPB後のILRMA推定信号(周波数領域)
%        >estimate(1).Y_ILRMA_SP        : マイク1のPB後の提案法推定信号(周波数領域) IxJxN
%        >estimate(1).T_signal_ILRMA    : マイク1のILRMA推定信号(時間領域)
%        >estimate(1).T_signal_ILRMA_SP : マイク1の提案法推定信号(時間領域) sig_length x N    
%       cost: convergence behavior of cost function in ILRMA (it+1 x 1)

ns = Param.ns;
M = Param.mic;
nb = Param.nb;
fftSize = STFT.fftSize;
shiftSize = STFT.shiftSize;
it = Param.it;
type = Param.type;

%音源のSTFT(IとJを確保)
switch STFT.type
    case 1
        [S_1,~] = STFT_1(musicSS(1,:)',fftSize, shiftSize,STFT.w_name);
    case 2
        S_1 = STFT_basic(musicSS(1,:)',STFT.analysis_window, STFT.shiftSize, STFT.frag);
end
I=size(S_1,1);
J=size(S_1,2);

% Short-time Fourier transform
X=zeros(I,J,M);
for m=1:M
    switch STFT.type
        case 1
           [X(:,:,m), window] = STFT_1(mix(:,m), fftSize, shiftSize,STFT.w_name);
        case 2
            X(:,:,m) = STFT_basic(mix(:,m),STFT.analysis_window,STFT.shiftSize,STFT.frag);
    end
end 
% Whitening
% Xwhite = whitening( X, ns ); % decorrelate input multichannel signal by applying principal component analysis
T = rand(I,nb,ns);
V = rand(nb,J,ns);

%%%%%%%%%%% ILRMA %%%%%%%%%%
%[Y, cost] = ILRMA( Xwhite, type, it, nb, drawConv, normalize );
[Y1,~,A1,~,~] = ILRMA1_initial(X,T,V,Param);
%[Y1,~,A1,~,~] = ILRMA1_kitamura_initial(X,T,V,Param);
% [Y,cost,~] = ILRMA_readable(Xwhite, type, it, nb, drawConv, normalize);   
%[Y, cost] = ILRMA_readable_student(Xwhite, type, it, nb, drawConv, normalize);
%%%%%%%%%% 提案法 %%%%%%%%%%
[Y_SP,A_SP] = ILRMA1_Sparse(X,T,V,Param,Sparse,STFT.fftSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Back projection (fixing scale ambiguity using reference microphone)
%Z = backProjection( Y, X(:,:,refMic) ); % scale-fixed estimated signal
%Z_2 = backProjection(Y, X(:,:,2));      %マイク２参照

%ILRMA
% for m = 1:M
%     estimate(m).Y_ILRMA = zeros(I,J,ns);
%     estimate(m).Y_ILRMA = backProjection(Y1,X(:,:,m));
% end
Z_PB = projectionBack(Y1,A1);
for m = 1:M
    estimate(m).Y_ILRMA = zeros(I,J,ns);
    for n = 1:ns
        estimate(m).Y_ILRMA(:,:,n) = Z_PB(:,:,n,m);
    end
end

%ILRMA+SP
Z_SP_PB = projectionBack(Y_SP,A_SP);
for m = 1:M
    estimate(m).Y_ILRMA_SP = zeros(I,J,ns);
    for n = 1:ns
        estimate(m).Y_ILRMA_SP(:,:,n) = Z_SP_PB(:,:,n,m);
    end
end

%%%%% ISTFT %%%%%
for m=1:M
    estimate(m).T_signal_ILRMA = zeros(Param.sig_length,Param.ns);
    estimate(m).T_signal_ILRMA_SP = zeros(Param.sig_length,Param.ns);
    
    switch STFT.type
        case 1
          estimate(m).T_signal_ILRMA = ISTFT(estimate(m).Y_ILRMA,shiftSize,window,size(mix,1));
          estimate(m).T_signal_ILRMA_SP = ISTFT(estimate(m).Y_ILRMA_SP,shiftSize,window,size(mix,1));
        case 2
          for n = 1:ns
          estimate(m).T_signal_ILRMA(:,n) = ISTFT_basic(estimate(m).Y_ILRMA(:,:,n),STFT.dual_window,STFT.shiftSize,Param.sig_length,STFT.frag);
          estimate(m).T_signal_ILRMA_SP(:,n) = ISTFT_basic(estimate(m).Y_ILRMA_SP(:,:,n),STFT.dual_window,STFT.shiftSize,Param.sig_length,STFT.frag);
          end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sep_p = permute(sep,[2,1,3]);      %マイク１における推定音 source x samples x channel
%sep_p2 = permute(sep_2,[2,1,3]);   %マイク２における推定音

% sig_resample = permute(sig_resample,[3,1,2]); %source*samples*channelの形に

 %[SDR1,SIR1,SAR1,~] = bss_eval_sources(sep_p,musicSS);
 %[SDR2,SIR2,SAR2,~] = bss_eval_sources(sep_p2,musicSS);
 
 %各マイクとの平均
 %SDR_a = (SDR1+SDR2)/2
 %SIR_a = (SIR1+SIR2)/2
 %SAR_a = (SAR1+SAR2)/2
