function SNR = wada_snr(D,SR,DEBUG)
% SNR = wada_snr(D,SR)
%    Blind estimation of speech signal SNR 
%    based on waveform amplitude distribution.
%    D is the waveform (or name of soundfile).  
%    SR is waveform's sample rate (not actually used).
%    SNR returns estimated SNR in dB.
%    From:
%      "Robust signal-to-noise ratio estimation based on waveform
%       amplitude distribution analysis", Chanwoo Kim & Richard M Stern
%      Interspeech 2008, 2598–2601.
%
% 2010-12-03 Dan Ellis dpwe@ee.columbia.edu for RATS

if nargin < 3; DEBUG = 0; end

if ischar(D)
  [D,SR] = wavread(D);
end

% Precomputed table of integral from paper (from Chanwoo's C++ implementation)
table = 'Alpha0.400000.txt';
[dbvals,a,b,Gvals] = textread(table,'%f %s %s %f');

aD = abs(D);
aD(aD < 1e-10) = 1e-10;
dVal1 = mean(aD);
dVal2 = mean(log(aD));
dEng = sum(D.^2);

if dVal1 == 0;	dVal1 = 1e-10;   end

dVal3 = log(dVal1) - dVal2;

if DEBUG; fprintf('E[|z|] = %f E[log|z|] =  %f\n', dVal1, dVal2); end
if DEBUG; fprintf('log(E[|z|]) - E[log(|z|)] = %f\n', dVal3); end

% Table interpolation
%dSNR = dbvals(min(find(dVal3<=Gvals)));
dSNRix = max(find(Gvals < dVal3));
if length(dSNRix) == 0
  dSNR = dbvals(1);
elseif dSNRix == length(dbvals)
  dSNR = dbvals(end);
else % can actually interpolate
  dSNR = dbvals(dSNRix) ...
         + (dVal3-Gvals(dSNRix))/(Gvals(dSNRix+1)-Gvals(dSNRix))...
           *(dbvals(dSNRix+1)-dbvals(dSNRix));
end


% Calculate SNR
dFactor = 10^(dSNR / 10);
dNoiseEng   = dEng / (1 + dFactor);
dSigEng = dEng * dFactor / (1 + dFactor);

SNR = 10 * log10(dSigEng / dNoiseEng);

if DEBUG; fprintf('The computed SNR value is %f.\n', dSNR); end
if DEBUG; fprintf('Signal energy in this block : %f.\n', dSigEng); end
if DEBUG; fprintf('Noise energy in this block : %f.\n', dNoiseEng); end
if DEBUG; fprintf('Computed SNR : %f.\n', SNR); end
