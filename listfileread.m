function [L,N] = listfileread(F)
% [L,N] = listfileread(F)   Read a list of per-line items
%    F is a file containing a list of items, one per line.
%    Return L as a cell array, with one item per line, N as the
%    number of items.
%    If F is not found, return empty L and N == -1 (instead of 0).
% 2006-08-06 dpwe@ee.columbia.edu for MIREX 06
% 2008-08-11 Dan Ellis dpwe@ee.columbia.edu  updated with "@include" syntax
% 2011-08-02 Modified to skip blank lines and lines starting with '#'

% disp(['listfileread: ',F]);

N = -1;
L = [];

inckey = '@include';
commentchar = '#';  % hash in first column means line is ignored

if fexist(F) == 1

  fid = fopen(F);

  nitems = 0;

  while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end

    if strncmp(tline, inckey, length(inckey))
      % Special syntax to embed another list
      pdir = fileparts(F);
      % read in the specified file
      subfname = tline((length(inckey)+2):end);
      sdir = fileparts(subfname);
      if length(sdir) == 0 | sdir(1) == '.'
        subfname = fullfile(pdir,subfname);
      end
      [LL,NN] = listfileread(subfname);
      for i = 1:NN
        nitems = nitems+1;
        L{nitems} = LL{i};
      end
    else
      if length(tline) > 0
        if tline(1) ~= commentchar
          nitems = nitems+1;
          L{nitems} = tline;
        end
      end
    end
  end
  fclose(fid);
  N = nitems;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = fexist(F)
%  E = fexist(F)  returns 1 if file F exists, else 0
% 2006-08-06 dpwe@ee.columbia.edu

x = dir(F);

E = length(x);
