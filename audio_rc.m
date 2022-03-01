function [mix,audio,Param] = audio_rc(audio,Param,STFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%音源の読み込みとインパルス応答の畳み込みを行います．
%STFTはここでは行いません．あくまで混合(観測)信号の作成．
%
%[Inputs]
%       audio: 各音源に関する情報を持つ構造体
%       Param: mainで定義した各パラメータの構造体
%        STFT: STFTに関するパラメータの構造体
%[Outputs]
%         mix:畳み込み後の混合信号 (mic x justsamples)
%       audio:音源情報を持つ構造体
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = Param.ns;
M = Param.mic;
fsResample = Param.fsResample;

%Input data for all sound source 
for i = 1:N
   [audio(i).sig, fs] = audioread(audio(i).name);
end

%Resample
for i = 1:N
   audio(i).resample = resample(audio(i).sig, fsResample, fs); 
end

%just_sample = fix((size(audio(1).resample,1)-fftSize)/shiftSize)*shiftSize+fftSize;
Param.sig_length = fix((fsResample*Param.seconds-STFT.fftSize)/STFT.shiftSize)*STFT.shiftSize+STFT.fftSize;

%Ajustment of resample length  
%Changing monoral signal
for i = 1:N
   audio(i).resample = audio(i).resample(1:Param.sig_length,1)';
end

mix = zeros(M,size(audio(1).resample,2));
for n = 1:N
audio(n).music = zeros(M,size(audio(1).resample,2));
end

%Convolution
for n = 1:N
  for m = 1:M
      music = conv(audio(n).resample,audio(n).impulse(m,:));
      audio(n).music(m,:) = music(1:length(audio(n).resample));
  end
end

%Make mixed signal
for m = 1:M
    for n = 1:N
      mix(m,:) = mix(m,:) + audio(n).music(m,:);
    end
end

end
%%%%%%%%%%EOF%%%%%%%%%%