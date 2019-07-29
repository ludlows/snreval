function S = nist_stnr(D,SR,FORCE_MATLAB)
% S = nist_stnr(D,SR,FORCE_MATLAB)
%    Calculate NIST STNR
%    Currently works by running external binary & parsing return
%    But will fall back on pure matlab approximation, or use it
%    whenever FORCE_MATLAB is set
% 2010-12-02 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; SR = 8000; end
if nargin < 3; FORCE_MATLAB = 0; end
% Accept a filename instead of a vector
if ischar(D)
  [D,SR] = audioread(D);
end

sndcat = findbinarch('sndcat');
stnr   = findbinarch('stnr');

if FORCE_MATLAB || exist(sndcat) == 0 || exist(stnr) == 0
  disp(['sndcat / stnr binaries not found/disabled, using Matlab ' ...
        'approximation']);
  S = nist_stnr_m(D,SR);
else
  % save the extracted region of the noisy file
  noisyfile = 'dn.wav';
%   wavwrite(D,SR,noisyfile);
  audiowrite(noisyfile, D, SR);
  % convert to sphere
  noisysph = 'dn.sph';
  mysystem([sndcat, ' -T NIST ',noisyfile,' -o ',noisysph]);
  % run nist stnr on it
  rs = tokenize(mysystem([stnr,' -c ',noisysph, ...
                      ' | grep SNR']));
  % extract the reported SNR
  S = str2num(rs{14});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = mysystem(cmd)
% Run system command; report error; strip all but last line
[s,w] = system(cmd);
if s ~= 0 
  error(['unable to execute ',cmd,' (',w,')']);
end
% Keep just final non-blank line
if ~isempty(w)
  while(w(end)==10); w = w(1:end-1); end
  w = w((1+max([0,findstr(w,10)])):end);
end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = findbinarch(n)
% Attempt to find an executable binary whose core is n.
% Return full path in p, or empty if not found.

% do we have the binaries?
arch = computer();
% fallback to 32 bit linux if 64 bit linux
if strcmp(arch,'GLNXA64') == 1
  arch = 'GLNX86';
end

p = findbin(n);
if isempty(p)
  p = findbin([n, '.', arch]);
end

%disp(['returning ',p]);

%%%%%%%%%%%%%%%%%%%%%%%
function p = findbin(n)
% Actually try "which", or local dir, to find bin named n
%disp(['trying ',n,' ...'])
[r, p] = system(['which ', n]);
if r ~= 0
  % Try here
  p = fullfile('.', n);
  %disp(['trying ',p,' ...'])
  r = (exist(p) == 0);
end
if r ~= 0
  p = '';
end

% strip any returns
p = p(double(p)>31);

%disp(['findbin: returning ',p]);

