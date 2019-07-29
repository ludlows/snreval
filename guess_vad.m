function T = guess_vad(D,SR,tsm)
% T = guess_vad(D,SR,tsm)
%     Estimate voice activity times from a waveform.
%     Return T as a set of [voice_start voice_end] rows, in
%     seconds.
%     tsm is median smoothing time, default 0.25 s.
% 2010-12-02 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; tsm = 0.25; end

if ischar(D)
  [D,SR] = wavread(D);
end

% Choose fft size; 256 for sr = 16000
nfft = 2^round(log(256 * SR/16000)/log(2));

% Make spectrogram
DC = 20*log10(abs(specgram(D,nfft)));
% Frame rate for spectrogram
fr = SR/(nfft/2);

% Take energy in voice region 100..1000 Hz
vxfmin = 100;
vxfmax = 1000;
leb = round(vxfmin*nfft/SR);
ueb = round(vxfmax*nfft/SR);
% voicing energy as largest bin in that range (?)
DCE = max(DC(leb:ueb,:));
% remove -Infs, or anything more than 50 dB below peak
lowthresh = max(DCE) - 50.0;
DCE(DCE < lowthresh) = lowthresh;

% then threshold as 1/3rd way between 10th and 90th percentile, in dB
%thr = percentile(DCE,0.1)+6.0; % twice the 10th percentile, in dB
thr = 0.33*percentile(DCE,0.1) + 0.67*percentile(DCE,0.9);

% Median smooth over 250ms
%tsm = 0.25;
%disp(['smoothing time = ',num2str(tsm)]);
% odd-length median filter window
medwin = 1+2*round(tsm/2*SR/(nfft/2));
% VAD estimate
vad = (DCE > thr);
% Blur out 1 bin each way
vad = max([vad;[0,vad(1:end-1)];[vad(2:end),0]]);
% Median filter
vad = medianf(vad, medwin);

% debug plot
debug =0;
if debug

  subplot(211)
  imgsc(DC)
  
  subplot(212)
  
  tt = [1:length(DCE)]*nfft/2/SR;
  ttt = tt([1 end]);
  plot(tt, DCE,'-b', ...
       ttt, [thr thr],'-r', ...
       ttt, percentile(DCE,0.1)*[1 1], ':g', ...
       ttt, percentile(DCE,0.9)*[1 1], ':r');
  title(['VAD threshold=',num2str(thr)])
   figure; % stop eval_snr overwriting this plot

end

% Actual times corresponding to sgram bins are hop/sr:hop/sr:end
tt = [1:length(DCE)]/fr;
%plot(tt,DCE,[tt(1) tt(end)],[thr thr],tt,10*vad, 'r');

% Convert to times
voiceon = find(([0,vad]==0)&([vad,0]==1));
voiceoff = find(([0,vad]==1)&([vad,0]==0));

T = zeros(length(voiceon),2);

for i = 1:length(voiceon)
  T (i,:) = tt([voiceon(i) min(voiceoff(i),length(DCE))]);
end
