function [Y,FS,NBITS,OPTS] = flacread(FILE,N,MONO,DOWNSAMP,DELAY)
% FLACREAD   Read FLAC audio file via use of external binaries.
%   Y = FLACREAD(FILE) reads a FLAC-encoded audio file into the
%     vector Y just like wavread reads a wav-encoded file (one channel 
%     per column).  Extension ".flac" is added if FILE has none.
%     Also accepts other formats of wavread, such as
%   Y = FLACREAD(FILE,N) to read just the first N sample frames (N
%     scalar), or the frames from N(1) to N(2) if N is a two-element vector.  
%   Y = FLACREAD(FILE,FMT) or Y = flacread(FILE,N,FMT) 
%     with FMT as 'native' returns int16 samples instead of doubles; 
%     FMT can be 'double' for default behavior (to exactly mirror the
%     syntax of wavread).
%
%   [Y,FS,NBITS,OPTS] = FLACREAD(FILE...) returns extra information:
%     FS is the sampling rate,  NBITS is the bit depth (always 16), 
%     OPTS.fmt is a format info string; OPTS has multiple other
%     fields, see WAVREAD.
%
%   SIZ = FLACREAD(FILE,'size') returns the size of the audio data contained
%     in the file in place of the actual audio data, returning the
%     2-element vector SIZ=[samples channels].
%
%   [Y...] = FLACREAD(FILE,N,MONO,DOWNSAMP) extends the
%     WAVREAD syntax to emulate special features of the
%     mpg123 engine:  MONO = 1 forces output to be mono (by
%     averaging stereo channels); DOWNSAMP = 2 or 4 downsamples by 
%     a factor of 2 or 4 (thus FS returns as 22050 or 11025
%     respectively for a 44 kHz flac file).  (A final DELAY argument
%     is also supported for full compatibility with mp3read, but
%     has no use).  In this case, N is interpreted in terms of the 
%     post-downsampling samples
%
%   Example:
%   To read an flac file as doubles at its original width and sampling rate:
%     [Y,FS] = flacread('piano.flac');
%   To read the first 1 second of the same file, downsampled by a
%   factor of 4, cast to mono, using the default filename
%   extension:
%     [Y,FS4] = flacread('piano', FS/4, 1, 4);
%
%   Note: requires external binary flac
%     http://labrosa.ee.columbia.edu/matlab/flacread.html
%
%   See also wavread, mp3read, m4aread.

% $Header: /Users/dpwe/matlab/columbiafns/m4aread/RCS/m4aread.m,v 1.1 2010/09/17 16:07:46 dpwe Exp dpwe $

% 2011-03-08 flacread created from m4aread

% find our baseline directory
path = fileparts(which('flacread'));

% %%%%% Directory for temporary file (if needed)
% % Try to read from environment, or use /tmp if it exists, or use CWD
tmpdir = getenv('TMPDIR');
if isempty(tmpdir) || exist(tmpdir,'file')==0
  tmpdir = '/tmp';
end
if exist(tmpdir,'file')==0
  tmpdir = '';
end
% ensure it exists
%if length(tmpdir) > 0 && exist(tmpdir,'file')==0
%  mkdir(tmpdir);
%end

%%%%%% Command to delete temporary file (if needed)
rmcmd = 'rm';

%%%%%% Location of the binaries - attempt to choose automatically
%%%%%% (or edit to be hard-coded for your installation)
ext = lower(computer);
if ispc
  ext = 'exe';
  rmcmd = 'del';
end
%flac = fullfile(path,['flac.',ext]);
%metaflac = fullfile(path,['metaflac.',ext]);
%flac = '/opt/local/bin/flac';
%metaflac = '/opt/local/bin/metaflac';
[r,flac] = system('which flac');
if r ~= 0; error(flac); end
[r,metaflac] = system('which metaflac');
if r ~= 0; error(metaflac); end
% strip trailing returns
flac = flac(1:end-1);
metaflac = metaflac(1:end-1);

%%%%% Process input arguments
if nargin < 2
  N = 0;
end

% Check for FMT spec (per wavread)
FMT = 'double';
if ischar(N)
  FMT = lower(N);
  N = 0;
end

if length(N) == 0
  N = [1 0];
elseif length(N) == 1
  % Specified N was upper limit
  N = [1 N];
end
if nargin < 3
  forcemono = 0;
else
  % Check for 3rd arg as FMT
  if ischar(MONO)
    FMT = lower(MONO);
    MONO = 0;
  end
  forcemono = (MONO ~= 0);
end
if nargin < 4
  downsamp = 1;
else
  downsamp = DOWNSAMP;
end
if downsamp ~= 1 && downsamp ~= 2 && downsamp ~= 4
  error('DOWNSAMP can only be 1, 2, or 4');
end

% process DELAY option (nargin 5) after we've read the SR

if strcmp(FMT,'native') == 0 && strcmp(FMT,'double') == 0 && ...
      strcmp(FMT,'size') == 0
  error(['FMT must be ''native'' or ''double'' (or ''size''), not ''',FMT,'''']);
end


%%%%%% Constants
NBITS=16;

%%%%% add extension if none (like wavread)
[path,file,ext] = fileparts(FILE);
if isempty(ext)
  FILE = [FILE, '.flac'];
end

  %%%%%% Probe file to find format, size, etc. 
  cmd=['"',metaflac,'" --show-total-samples --show-sample-rate ', ...
       '--show-channels --show-bps "',FILE,'"'];
  [s,w] = system(cmd);
  if s ~= 0 
    error(['unable to execute ',cmd,' (',w(1:end-1),')']);
  end
  xx = find(w==10);
  nframes = str2num(w(1:xx(1)-1));
  SR = str2num(w((xx(1)+1):(xx(2)-1)));
  nchans = str2num(w((xx(2)+1):(xx(3)-1)));
  NBITS = str2num(w((xx(3)+1):(xx(4)-1)));
  
  smpsperfrm = nchans;

  % fields from wavread's OPTS
  OPTS.fmt.nAvgBytesPerSec = 0; % bitrate/8;
  OPTS.fmt.nSamplesPerSec = SR;
  OPTS.fmt.nChannels = nchans;
  OPTS.fmt.nBlockAlign = 0; %smpspfrm/SR*bitrate/8;
  OPTS.fmt.nBitsPerSample = NBITS;
  
% process or set delay
if nargin < 5
  delay = 0;
else
  delay = DELAY;
end

% Size-reading version
if strcmp(FMT,'size') == 1
   Y = [floor(nframes/downsamp), nchans];
   FS = SR;
else

  % Temporary file to use
  %tmpfile = fullfile(tmpdir, ['tmp',num2str(round(1000*rand(1))),'.wav']);
  [s,upid] = system('echo $$');
  % remove final CR
  upid = upid(1:end-1);
  tmpfile = fullfile(tmpdir, ['flactmp',upid,'.wav']);

  durstr = '';
  sttfrm = N(1)-1;
  if sttfrm > 0
    durstr = ['--skip=',num2str(sttfrm*downsamp)];
  end
  endfrm = 0;
  if length(N) > 1
    endfrm = N(2);
    if endfrm > 0
      durstr = [durstr,' --until=',num2str(endfrm*downsamp),];
    end
  end
  

  % Run the decode
  cmd=['"',flac,'" -o "', tmpfile,'" ',durstr,' -d "',FILE,'"'];
  %w = 
  mysystem(cmd);

  % Load the data
  [Y,SR] = wavread_downsamp(tmpfile, [], forcemono, downsamp);

  % Delete tmp file
  mysystem([rmcmd,' "', tmpfile,'"']);
  
  % debug
%  disp(['sttfrm=',num2str(sttfrm),' endfrm=',num2str(endfrm),' skipx=',num2str(skipx),' delay=',num2str(delay),' len=',num2str(length(Y))]);
  
  % Convert to int if format = 'native'
  if strcmp(FMT,'native')
    Y = int16((2^15)*Y);
  end

  FS = SR;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = mysystem(cmd)
% Run system command; report error; strip all but last line
[s,w] = system(cmd);
if s ~= 0 
  error(['unable to execute ',cmd,' (',w(1:end-1),')']);
end
% Keep just final line
w = w((1+max([0,findstr(w,10)])):end);
% Debug
%disp([cmd,' -> ','*',w,'*']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = tokenize(s,t)
% Break space-separated string into cell array of strings.
% Optional second arg gives alternate separator (default ' ')
% 2004-09-18 dpwe@ee.columbia.edu
if nargin < 2;  t = ' '; end
a = [];
p = 1;
n = 1;
l = length(s);
nss = findstr([s(p:end),t],t);
for ns = nss
  % Skip initial spaces (separators)
  if ns == p
    p = p+1;
  else
    if p <= l
      a{n} = s(p:(ns-1));
      n = n+1;
      p = ns+1;
    end
  end
end
    
