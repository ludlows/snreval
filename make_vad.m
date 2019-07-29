function make_vad(F,V,G)
% make_vad(F,V,G)
%     Read clean audio file F, run guess_vad on it, write it out to
%     V (defaults to F stem + '-vad.txt');
%     G = 1 means plot a graph too.
% 2011-03-29 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2;  V = ''; end

if nargin < 3;  G = 0; end

if length(V) == 0
  [p,n,e] = fileparts(F);
  V = fullfile(p, [n,'-vad.txt']);
end

[d,sr] = audioread(F);
vadtimes = guess_vad(d,sr);

if G
  plot([1:length(d)]/sr, d);
  hold on; plot(vadtimes', 0*vadtimes'+.2, '-r'); hold off
  hold on; plot(vadtimes', 0*vadtimes'-.2, '-r'); hold off
  hold on; plot([1;1]*[vadtimes(:,1)',vadtimes(:,2)'], [0.2;-0.2]*ones(1,size(vadtimes,1)*2), '-r'); hold off
end

write_vad_file(vadtimes, V);
disp(['Wrote ',num2str(size(vadtimes,1)),' voicing events to ',V]);