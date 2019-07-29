function write_vad_file(T,F)
% write_vad_file(T,F)
%   Write a voice activity file F from the times in T, each row is 
%   [voicestarttime voiceendtime] in secs.
% 2010-12-02 Dan Ellis dpwe@ee.columbia.edu

f = fopen(F,'w');
fprintf(f,'%f %f 1\n',T');
fclose(f);
