function [SNRstnro,SNRwada,SNRvad,SAR,pesqmos,targdelay] = snreval(NOISY, varargin)
% [SNRstnr,SNRwada,SNRvad,SAR,pesqmos,targdelay] = snreval(noisyfile, varargin)
%    Run various measures on a file to estimate SNR.
%  NOISY is the name of the noisy file.
%  '-vad', VAD gives the name of the hand-marked voice activity file
%  '-clean', CLEAN gives the name of a corresponding clean-speech file
%  '-start', TS, '-end', TE specify start and end times within
%    the NOISY file for processing a subrange (also applied to
%    CLEAN file, if any). 
%  '-guessvad', 1 - try to guess the VAD from CLEAN (or NOISY if
%    no CLEAN).

%  '-disp', 0 - don't do any graphics.
%  '-listin', 1 - treat NOISY as a text file listing the actual
%    noisy files to process
%  '-listout', 1, - write output values in columns instead of text report
%  '-vaddir', DIR - for each noisy file, read a corresponding file
%    from this directory containing the VAD information
%  '-vadlist', FILE - file containing list of VAD files
%  '-cleandir', DIR - for each noisy file, read a corresponding file
%    from this directory containing the clean signal
%  '-cleanlist', FILE - file containing list of CLEAN files
%  '-ldclabels', 1 - treat VAD file as 8-column LDC format
%  '-samplerate', SR - resample data to this SR before processing
%  '-checkfshift', 1 - search for SSB frequency shift too
%  '-hpf', FREQ - high-pass input signals at this frequency (for LF noise)
%  '-preemp', MAG - pre-emphasis zero magnitude (e.g. 0.97, default 0 = none)
%  '-my_stnr', 1 - force pure-matlab version of STNR (not nist binary)
%
%    For RATS.
% 2010-12-14 Dan Ellis dpwe@ee.columbia.edu

VERSION = 0;
DATE = 20140701;

% Parse out the optional arguments
[VAD,CLEAN,TS,TE,GUESS,DISP,LISTIN,LISTOUT, ...
 VADDIR,VADLIST,CLEANDIR,CLEANLIST,LDCLABELS,SAMPLERATE, ...
 CHECKFSHIFT, HPF, PREEMPH, MY_STNR, ...
 XTRA] = ...
    process_options(varargin, '-vad', '', '-clean', '', ...
                    '-start', 0, '-end', 0', '-guessvad', 0, ...
                    '-disp', 1, ...
                    '-listin', 0, '-listout', 0, ...
                    '-vaddir', '', '-vadlist', '', ...
                    '-cleandir', '', '-cleanlist', '', ...
                    '-ldclabels', 0, '-samplerate', 0, ...
                    '-checkfshift', 0, '-hpf', 0, '-preemph', 0, ...
                    '-my_stnr', 0);

HELP = 0;
if length(XTRA) > 0
  HELP = length(strmatch('-help',XTRA,'exact')) > 0;
  if ~HELP
    disp(['Unrecognized options:',sprintf(' %s',XTRA{1:end})]);
    HELP = 1;
  end
end

if strcmp(NOISY,'-help') == 1; HELP = 1; end
  
if HELP
  disp(['snreval v',num2str(VERSION),' of ',num2str(DATE)]);
  help('snreval');
  return
end

if LISTOUT || ~DISP
  QUIET = 1;
else
  QUIET = 0;
end


% Text header
disp(['#============= SNREVAL v',num2str(VERSION),' (',num2str(DATE),') ===']);
disp(['# args: ', join(varargin)]);
% If LISTOUT, print out column headers
if LISTOUT
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
               '#TARGETFILE', 't_start', 't_end', 't_delay', ...
               'STNR', 'WADA', 'SNRvad', 'SAR', 'PESQ');
end


do_guess_vad = GUESS;

% order of high-pass filters
HPFORD = 8;

SNRvad = 0;
SAD = 0;
pesqmos = 0;
targdelay = 0;

% default output, so doesn't crash if input list is empty
SNRstnr = [];

if LISTIN
  noisylist = listfileread(NOISY);
else
  noisylist = {NOISY};
end

if length(VADLIST) > 0
  vadlist = listfileread(VADLIST);
else
  vadlist = '';
end
if length(CLEANLIST) > 0
  cleanlist = listfileread(CLEANLIST);
else
  cleanlist = '';
end

nnoisy = length(noisylist);

for nf = 1:nnoisy
  
  % initialize values needed for list output
  targdelay = 0;
  SNRstnr = -999;
  SNRwada = -999;
  SNRvad = -999;
  SAR = -999;
  pesqmos = -999;

  
  NOISY = noisylist{nf};
  [p,n,e] = fileparts(NOISY);
  
  % maybe construct VAD and CLEAN files
  if length(VADDIR) > 0
    VAD = fullfile(VADDIR, [n,'.txt']);
  end
  if length(CLEANDIR) > 0
    CLEAN = fullfile(CLEANDIR, [n,e]);
  end
  if length(vadlist) > 0
    VAD = vadlist{nf};
  end
  if length(cleanlist) > 0
    CLEAN = cleanlist{nf};
  end

  % Read in the specified VAD file, if any
  if length(VAD) > 0
    [vadtimes,vnatimes] = read_vad_file(VAD, LDCLABELS);
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
  DUR = max(0,TE-TS);
%   [dn,sr] = audioread(NOISY,SAMPLERATE,1,TS,DUR);
  [dn,sr] = audioread(NOISY);
  if HPF > 0
    [bh,ah] = butter(HPFORD, HPF/(sr/2), 'high');
    dn = filter(bh, ah, dn);
  end
  if PREEMPH ~= 0
    dn = filter([1 -PREEMPH], 1, dn);
  end
  
  dlen = length(dn);

  [pp,nm,ee] = fileparts(NOISY);
  nm(nm=='_') = ' ';

  if DISP
    % Reset all the subplots
    subplot(111)
    plot(1)
  end
  
  if have_clean
%     [dc,src] = audioread(CLEAN,sr,1,TS,DUR);
    [dc,src] = audioread(CLEAN);
    assert(src==sr);
    if HPF > 0
      [bh,ah] = butter(HPFORD, HPF/(sr/2), 'high');
      dc = filter(bh, ah, dc);
    end
    if PREEMPH ~= 0
      dc = filter([1 -PREEMPH], 1, dc);
    end
    dlen = min(dlen,length(dc));
  end

  if TE == 0 || TE > TS+dlen/sr
    thisTE = TS+dlen/sr;
  else
    thisTE = TE;
  end

  if thisTE > 0 || thisTE < TS+dlen/sr
    % truncate files as indicated
    dn = dn(1:round((thisTE-TS)*sr));
    dlen = length(dn);
    if have_clean
      dc = dc(1:round((thisTE-TS)*sr));
      dlen = min(dlen,length(dc));
    end
  end

  if have_vad == 0
    if do_guess_vad > 0
      if have_clean
        if ~QUIET
          disp(['Guessing VAD from clean file ', CLEAN,' ...']);
        end
        vadsource = dc;
      else
        if ~QUIET
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
      vnatimes = [vadtimes(1:end-1,2),vadtimes(2:end,1)];
      have_vad = 1;
    else
      if ~QUIET
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
    vv = zeros(1,length(dn))==1;
    % and also sgram frames (separate masks for active and inactive)
    vxf = zeros(1,nsgframes);
    nvf = zeros(1,nsgframes);
    lastvoff = 1+0;
    for i = 1:size(vadtimes,1)
      von = vadtimes(i,1);
      vof = vadtimes(i,2);
      % skip any entries entirely outside our time range
      if vof > TS && von < thisTE
        % map into our time range & clip
        von = max(0,von - TS);
        vof = min(thisTE-TS,vof - TS);
        vv((1+round(von*sr)):(round(vof*sr))) = (1==1);
        % only mark as voiced frames that are completely within voiced region
        vonf = von*fr;
        voff = vof*fr;
        vxf(ceil(1+vonf):floor(voff)) = 1;
      end
    end 
    for i = 1:size(vnatimes,1)
      non = vnatimes(i,1);
      nof = vnatimes(i,2);
      % skip any entries entirely outside our time range
      if nof > TS && non < thisTE
        % map into our time range & clip
        non = max(0,non - TS);
        nof = min(thisTE-TS,nof - TS);
        % only mark as unvoiced frames that are completely in unvoiced region
        nonf = non*fr;
        noff = nof*fr;
        nvf(ceil(1+nonf):floor(noff)) = 1;
      end
    end 
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
  SNRstnr = nist_stnr(dn,sr,MY_STNR);
  msgstnr = ['NIST STNR = ',mynum2str(SNRstnr),' dB'];

  % Calculate the WADA SNR
  SNRwada = wada_snr(dn,sr);
  msgwada = ['WADA SNR  = ',mynum2str(SNRwada),' dB'];

  if DISP
    % Plot waveform and VAD
    subplot(411)
    plot([1:length(dn)]/sr, dn);
    if have_vad
      hold on; plot([0:length(vv)-1]/sr,0.2*vv,'-g'); hold off
      hold on; plot([0:length(vv)-1]/sr,-0.2*vv,'-g'); hold off
    end
    title(['File: ',nm,' - ',mynum2str(TS),'-',mynum2str(thisTE),' s']);
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
      % Overplot vad on spectrogram, first notactive, then active
      fmin = 100;
      fmax = 0.9*sr/2;
      for i = 1:size(vnatimes,1)
        non = vnatimes(i,1)-TS;
        nof = vnatimes(i,2)-TS;
        hold on; 
        plot([non nof nof non non],...
             [fmin fmin fmax fmax fmin],'-r'); 
        hold off
      end
      for i = 1:size(vadtimes,1)
        von = vadtimes(i,1)-TS;
        vof = vadtimes(i,2)-TS;
        hold on; 
        plot([von vof vof von von],...
             [fmin fmin fmax fmax fmin],'-g'); 
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
      title(['SNRvad    = ',mynum2str(SNRvad),' dB']);
    end

    msgsnrv = ['   SNRvad = ',mynum2str(SNRvad),' dB'];
  end

  if have_clean

    % Calculate PESQ
    % save out the clean file (noisy already written, but do it again
    % to be safe)
    noisyfile = 'dn.wav';
%     wavwrite(dn,sr,noisyfile);
    audiowrite(noisyfile, dn, sr);
    cleanfile = 'dc.wav';
%     wavwrite(dc,sr,cleanfile);
    audiowrite(cleanfile, dc, sr);
    % run the precompiled PESQ routine
    pesqlqo = pesq_wrap(cleanfile,noisyfile);
    msgpesqlqo = [' PESQ NARROW MOS-LQO = ',mynum2str(pesqlqo)];
    msgpesq = [' PESQ NARROW MOS     = ',mynum2str(mos2pesq(pesqlqo))];
  
    % Linear decomposition into target + artifacts
    % delay noise rel to clean a little more to permit warmup of filter
%    warmup = round(0.010*sr);
%    % we filt a filter of up to 50ms IR, and do it in 2s blocks
%    [s_targ, e_artif] = decomp_lin_win(dn(1:end-warmup),dc(warmup+1:end),...
%                                       round(0.05*sr),round(1.0*sr));

    if CHECKFSHIFT
      % checking for frequency shift is slow, so only do it if requested
      [e_artif, s_targ, dfilt, SNRin, targdelay, fshift] = ...
          find_in_mix_tf(dn, dc, sr);
      if ~QUIET
        disp(['Mix freq shift= ',sprintf('%.1f',fshift),' Hz']);
      end
    else
      [e_artif, s_targ, dfilt, SNRin, targdelay] = ...
          find_in_mix(dn, dc, sr);
      fshift = 0;
    end

    
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
  

  if LISTOUT
    fprintf('%s\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
                 NOISY, TS, thisTE, targdelay, SNRstnr, SNRwada, SNRvad, ...
                 SAR, pesqmos);
  else
    disp('==========================================')
    disp(['Target File: ',NOISY]);
    disp([' time range: ',mynum2str(TS),'-',mynum2str(thisTE),' s'])
    if have_clean 
      disp(['   Ref File: ',CLEAN]);
      disp([' Targ delay: ',sprintf('%.3f',targdelay),' s']);
    end
%    if HPF > 0
%      disp(['  high-pass: ',sprintf('%.0f Hz order %d Butterworth', ...
%                                    HPF,HPFORD)]);
%    end
    disp(msgstnr)
    disp(msgwada)
    if have_vad
      disp(msgsnrv)
    end
    if have_clean 
      disp(msgsar);
      disp(msgpesqlqo);
      disp(msgpesq);
      
    end
    disp('==========================================')
  end

end


  
%if DISP
%  disp('Press any key to continue...');
%  pause;
%end

% only set the output if there are some vars assigned
if nargout > 0
  SNRstnro = SNRstnr;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = mynum2str(n)
s = sprintf('%.5f',n);
end


function s = join(array)
s = '';
    for i = 1:length(array)
      val = array{i};
      if isstr(i)
        s = [s, val, ' '];
      else
        s = [s, num2str(val), ' '];
      end
    end
end

end