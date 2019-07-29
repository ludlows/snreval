function mos = pesq_wrap(clean, noisy)
% mos = pesq_wrap(clean, noisy)
%   Wrapper to execute Loizou's PESQ implementation, either from M
%   files or from pfiles
% 2011-02-24 Dan Ellis dpwe@ee.columbia.edu
[clean_wav, fs] = audioread(clean);
[noisy_wav, ~] =  audioread(noisy);
mos = pesq_mex(clean_wav,noisy_wav,fs, 'narrowband'); % wideband narrowband

end
