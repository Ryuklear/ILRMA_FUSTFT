function [audio] = reading_impulse(audio,impulse_name,Param,fftSize)
%RWCP音声音響データベース　E2Aインパルス応答データ
%データ　48kHz 32bit(float) rawファイル(バイナリデータ)
%E2Aは残響時間0.3sの残響可変室で収録されたデータ
%音源数×マイク数のインパルス応答を使用
Fs = 48000;

for n = 1:Param.ns
for m = 1:Param.mic
   fileID = fopen(impulse_name((n-1)*Param.mic+m,:));    
   
   h = fread(fileID,'float');
   %h = fread(fileID,'float');
   h = resample(h,Param.fsResample,Fs);
   fclose(fileID);
   audio(n).impulse(m,:) = h(1:fftSize);
end
end

%読み込んだデータは構造体audio.impulseとして返される

end
%%%%%EOF%%%%%