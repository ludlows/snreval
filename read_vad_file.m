function [TA,TI] = read_vad_file(F,FMT)
% [TA,TI] = read_vad_file(F,FMT)
%     Read a voice activity annotation file named F.
%     Return TA as a set of rows [voice_start voice_end] in seconds.
%     Return TI as a set of rows [novoice_start novoice_end]
%     to allow gaps that are neither.
%     FMT is 0 for "start end transcript" lines (default)
%     or 1 for 8-column LDC format (S/NS/NT)
% 2010-12-02 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2;  FMT = 0; end

if FMT ==0
  
  % Files written by Audacity as t_start t_end label
  [ss, ee, lab] = textread(F,'%f %f %s');

  %T = [ts(la>0),te(la>0)];
  TA = [ss,ee];
  TI = [ee(1:end-1),ss(2:end)];
  
else
  
  % FMT 1 - LDC 8 column
  % LDC2011E32 tab-separated columns
  % 1:filename 2:channel 3:starttime 4:endtime 
  % 5:label(S/NS/NT) 6:labelsource 7:language 8:langsource
  [fn,ch,ss,ee,lab,lsrc,lid,lsrc] = textread(F, ...
                             '%s\t%s\t%f\t%f\t%s\t%s\t%s\t%s', ...
                             'delimiter', '\t');
  labS = strmatch('S',lab,'exact');
  labNS = strmatch('NS',lab,'exact');

  TA = [ss(labS),ee(labS)];
  TI = [ss(labNS),ee(labNS)];
  
end
