function [SNRstnr,SNRwada,SNRvad,SAR,pesqmos,targdelay] = eval_snr(NOISY, VAD, CLEAN, TS, TE, GUESS, DISP)
% [SNRstnr,SNRwada,SNRvad,SAR,pesqmos,targdelay] = eval_snr(NOISY, VAD, CLEAN, TS, TE, GUESS, DISP)
%    Run various measures on a file to estimate SNR.
%    NOISY is the name of the noisy file.
%    VAD is the name of the hand-marked voice activity file
%    CLEAN is the name of a corresponding clean-speech file
%    TS, TE specify start and end times within the NOISY file for
%    processing a subrange (also applied to CLEAN file, if any).
%    If GUESS is set, try to guess the VAD from CLEAN (or NOISY if
%    no CLEAN).
%    If DISP is set to zero, don't do any graphics.
%    For RATS.
% 2010-11-30 Dan Ellis dpwe@ee.columbia.edu

VERSION = 0.3;
DATE = 20110802;

if nargin < 2;  VAD = ''; end
if nargin < 3;  CLEAN = ''; end
if nargin < 4;  TS = 0; end
if nargin < 5;  TE = 0; end
if nargin < 6;  GUESS = 0; end
if nargin < 7;  DISP = 1; end

% allow string inputs (for launching from command line)
if ischar(TS); TS = str2num(TS); end
if ischar(TE); TE = str2num(TE); end
if ischar(GUESS); GUESS = str2num(GUESS); end

do_guess_vad = GUESS;

% Allow both VAD and CLEAN files to be non-filenames as 
% an alias for leaving them empty -- because it's hard 
% to pass empty strings as arguments through the MCR
% shell script interface
if length(VAD) > 0 && exist(VAD,'file') == 0
  disp(['VAD file "',VAD,'" not found, ignoring...']);
  VAD = '';
end
if length(CLEAN) > 0 && exist(CLEAN,'file') == 0
  disp(['CLEAN file "',CLEAN,'" not found, ignoring...']);
  VAD = '';
end

SNRvad = 0;
SAD = 0;
pesqmos = 0;
targdelay = 0;

% Read in the specified VAD file, if any
if length(VAD) > 0
  vadtimes = read_vad_file(VAD);
  have_vad = 1;
else
  have_vad = 0;
end

if length(CLEAN) > 0
  have_clean = 1;
else
  have_clean = 0;
end

% Read in sound file
[dn,sr] = audioread(NOISY,[],1);
dlen = length(dn);

[pp,nm,ee] = fileparts(NOISY);
nm(nm=='_') = ' ';

if DISP
  % Reset all the subplots
  subplot(111)
  plot(1)
end
  
if have_clean
  [dc,src] = audioread(CLEAN,[],1);
  if src ~= sr
    dc = resample(dc,sr,src);
  end
  dlen = min(dlen,length(dc));
end

if TE == 0 || TE > dlen/sr
  TE = dlen/sr;
end

if TE > 0 || TE < dlen/sr
  % truncate files as indicated
  dn = dn((round(TS*sr)+1):round(TE*sr));
  dlen = length(dn);
  if have_clean
    dc = dc((round(TS*sr)+1):round(TE*sr));
    dlen = min(dlen,length(dc));
  end
end

if have_clean
  % find offset between clean and noisy
  xc = xcorr(dn, dc, sr/2);
  xcorig = find(abs(xc)==max(abs(xc)))-(sr/2+1);
  targdelay = xcorig/sr;
  
  % trim files to be aligned
  if xcorig > 0
    dn = dn(xcorig:end);
    TS = TS + xcorig/sr;
  elseif xcorig < 0
    dc = dc((-xcorig):end);
  end
  % make files same length
  dlen = min(length(dn),length(dc));
  dn = dn(1:dlen);
  dc = dc(1:dlen);
  TE = TS + dlen/sr;

  if DISP
    disp(['Optimal alignment: noisy is delayed by ',num2str(xcorig/sr),' s']);
  end
end

if have_vad == 0
  if do_guess_vad > 0
    if have_clean
      if DISP
        disp(['Guessing VAD from clean file ', CLEAN,' ...']);
      end
      vadsource = dc;
    else
      if DISP
        disp(['Guessing VAD from noisy file ', NOISY,' ...']);
      end
      vadsource = dn;
    end
    if do_guess_vad ~= 1.0
      vad_tsm = do_guess_vad;
    else
      vad_tsm = 0.25;
    end
    vadtimes = TS + guess_vad(vadsource, sr, vad_tsm);
    have_vad = 1;
  else
    if DISP
      disp('No VAD, but guessing is not selected');
    end
  end
end
 
% Choose fft size; 1024 for sr = 16000
nfft = 2^round(log(1024 * sr/16000)/log(2));
nhop = nfft/2;
nsgframes = 1+ floor((dlen-nfft)/nhop);
fr = sr/nhop;

if have_vad
  % Convert VAD times to indices
  % sample level
  vv = zeros(1,length(dn));
  % and also sgram frames (separate masks for active and inactive)
  vxf = zeros(1,nsgframes);
  nvf = zeros(1,nsgframes);
  lastvoff = 1+0;
  for i = 1:size(vadtimes,1)
    von = vadtimes(i,1);
    vof = vadtimes(i,2);
    % skip any entries entirely outside our time range
    if vof > TS && von < TE
      % map into our time range & clip
      von = max(0,von - TS);
      vof = min(TE-TS,vof - TS);
      vv((1+round(von*sr)):(round(vof*sr))) = 1;
      % only mark as voiced frames that are completely within voiced region
      vonf = von*fr;
      voff = vof*fr;
      vxf(ceil(1+vonf):floor(voff)) = 1;
      nvf(ceil(1+lastvoff):floor(vonf)) = 1;
      lastvoff = voff;
    end
  end 
  % convert vv to logical
  vv = (vv==1);
  % convert frame masks to indices
  vxf = find(vxf(1:nsgframes));
  nvf = find(nvf(1:nsgframes));
else
  % no VAD
  vv = ones(1,length(dn))==1;  % to make a logical
  vxf = 1:nsgframes;
  nvf = [];
end



% Calculate NIST STNR
% save the extracted region of the noisy file
SNRstnr = nist_stnr(dn,sr);
msgstnr = ['NIST STNR = ',mynum2str(SNRstnr),' dB'];

% Calculate the WADA SNR
SNRwada = wada_snr(dn,sr);
msgwada = ['WADA SNR  = ',mynum2str(SNRwada),' dB'];

if DISP
  % Plot waveform and VAD
  subplot(411)
  plot([1:length(dn)]/sr, dn);
  if have_vad
    hold on; plot([0:length(vv)-1]/sr,0.2*vv,'r'); hold off
    hold on; plot([0:length(vv)-1]/sr,-0.2*vv,'r'); hold off
  end
  title(['File: ',nm,' - ',mynum2str(TS),'-',mynum2str(TE),' s']);
end
  
% Plot signal spectrogram
DN = 20*log10(abs(specgram(dn,nfft)));
% collapse to 0..80 dB range (from top)
maxdn = max(DN(:));
DN = max(0,80+DN-max(DN(:)));

if DISP
  % Plot spectrogram
  subplot(412)
  frqs = [0:size(DN,1)-1]*sr/nfft;
  tts = [1:size(DN,2)]/fr;
  imgsc(tts,frqs,DN);
  caxis([0 80])
  title(['Noisy signal']);

  if have_vad
    % Overplot vad on spectrogram
    fmin = 100;
    fmax = 0.9*sr/2;
    for i = 1:size(vadtimes,1)
      von = vadtimes(i,1)-TS;
      vof = vadtimes(i,2)-TS;
      hold on; 
      plot([von vof vof von von],...
           [fmin fmin fmax fmax fmin],'-r'); 
      hold off
    end
  end
end

% Calculate subband energy histogram in each row, in 1 dB bins
dbmin = 0.0;
dbmax = 80.0;
dbbin = 1.0;
dbrange = [dbmin:dbbin:dbmax];
DNh = zeros(size(DN,1),length(dbrange));
for i = 1:size(DN,1);
  DNh(i,:) = hist(DN(i,:),dbrange);
end

if DISP
  subplot(425)
  imgsc(frqs,dbrange,DNh');
  maxdnh = max(max(DNh(:,2:end)));
  caxis([0 maxdnh/2])
  title(['Subband histograms - ', msgstnr,'; ', msgwada])
end

% .. and plot average of energy in voiced and unvoiced regions from VAD
if have_vad
  DNV = DN(:,vxf);
  DNNV = DN(:,nvf);
  Esn = (mean(sum(2*(10.^(DNV/10)))/nfft)/(nfft/2));
  En = (mean(sum(2*(10.^(DNNV/10)))/nfft)/(nfft/2));
  SNRvad = 10*log10(abs(Esn-En)/En);

  if DISP
    % Maybe plot actual spectrogram segments?
    plot_DNV = (have_clean==0);
    if plot_DNV
      subplot(4,2,7)
      %specgram(dn(vv==1),nfft,sr);
      imgsc(DNV);
      caxis([0 80])
      title('Spectrogram - voiced segments');
      subplot(4,2,8)
      %specgram(dn(nn==1),nfft,sr);
      imgsc(DNNV);
      caxis([0 80])
      title('Spectrogram - unvoiced segments');
    end
    % Plot averages in any case
    subplot(426)
    plot(frqs,mean(DNV,2),frqs,mean(DNNV,2),'-r')
    axis([min(frqs) max(frqs) 0 80])
    grid
    legend('Vx active','Noise')
    title(['SNRvad = ',mynum2str(SNRvad),' dB']);
  end

  msgsnrv = ['   SNRvad = ',mynum2str(SNRvad),' dB'];
end

if have_clean

  % Calculate PESQ
  % save out the clean file (noisy already written, but do it again
  % to be safe)
  noisyfile = 'dn.wav';
  wavwrite(dn,sr,noisyfile);
  cleanfile = 'dc.wav';
  wavwrite(dc,sr,cleanfile);
  % run the precompiled PESQ routine
  pesqmos = pesq_wrap(cleanfile,noisyfile);
  msgpesq = [' PESQ MOS = ',mynum2str(pesqmos)];
  
  % Linear decomposition into target + artifacts
  % delay noise rel to clean a little more to permit warmup of filter
  warmup = round(0.010*sr);
  % we filt a filter of up to 50ms IR, and do it in 2s blocks
  [s_targ, e_artif] = decomp_lin_win(dn(1:end-warmup),dc(warmup+1:end),...
                                     round(0.05*sr),round(1.0*sr));

  %SAR = 20*log10(std(s_targ)/std(e_artif));
  % Calculate SAR only over VAD active regions
  % Make sure files are long enough
  s_targ(length(vv)) = 0;
  e_artif(length(vv)) = 0;
  SAR = 20*log10(std(s_targ(vv))/std(e_artif(vv)));

  if DISP
    subplot(414)
    specgram(e_artif,nfft,sr);
    cax = caxis;
    %caxis([-80 0]+cax(2));
    caxis([-80 0]+maxdn);

    title(['Artifact residual - SAR = ',mynum2str(SAR),' dB, ', ...
           msgpesq])
  end

  msgsar = ['      SAR = ',mynum2str(SAR),' dB'];

end

if DISP
  colormap(1-gray);
end
  
%text(0,-0.5*(sr/2),msgstnr);
%text(dur/4,-0.5*(sr/2),msgsnrv);
%text(2*dur/4,-0.5*(sr/2),msgsar);
%text(3*dur/4,-0.5*(sr/2),msgpesq);


% Text report
disp(['============== SNREVAL v',num2str(VERSION),' (',num2str(DATE),') ===']);
disp(['Target File: ',NOISY]);
disp([' time range: ',mynum2str(TS),'-',mynum2str(TE),' s'])
if have_clean; 
  disp(['   Ref File: ',CLEAN]);
  disp([' Targ delay: ',sprintf('%.3f',targdelay),' s']);
end
disp(msgstnr)
disp(msgwada)
if have_vad
  disp(msgsnrv)
end
if have_clean; 
  disp(msgsar)
  disp(msgpesq)
end
disp('==========================================')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = mynum2str(n)
s = sprintf('%.1f',n);

