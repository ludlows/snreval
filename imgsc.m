function h = imgsc(a1, a2, a3)
% imgsc(X,Y,D)             Display a 2-D array of values as a pseudocolor plot.
%	Like imagesc, except vertical axis is flipped to make it more like 
%	pcolor(D).  X and Y specify the axis labelling, and can be omitted.
%	Handling of Y labels is funky and doesn't really work.
% dpwe 1994apr15.  Uses built-in 'imagesc'

v = version;
if(v([1 2 3])=='4.1')
  OLDVERS = 1;
else
  OLDVERS = 0;
end
if(v([1 2])=='5.' | v([1 2])=='6.' | v([1 2]) == '7.')    
  MAT5 = 1;
else
  MAT5 = 0;
end


if(nargin>=3)  % 3-arg version
  d = a3;
  s = size(d);
  t = size(a2);
  if(MAT5)
	  h = imagesc(a1, a2, d);
  else
	if(OLDVERS)
      h = imagesc(a1, a2([(t(1)*t(2)):-1:1]), d(s(1):-1:1,:));
  	else
% Apparently, because the second arg goes from smaller to larger, *don't* 
% need to flip matrix (4.2c-ism??)
      h = imagesc(a1, a2([(t(1)*t(2)):-1:1]), d);
	end
  end
else 
  d = a1;
  s = size(d);
  xx = 1:s(2);
  if(MAT5)
	yy = 1:s(1);
	h = imagesc(xx,yy,d);
  else
    yy = s(1):-1:1;
    if(OLDVERS)
      h = imagesc(d(s(1):-1:1,:));
    else
      h = imagesc(xx,yy,d);
    end
  end
end

% Matlab 5 way to get (0,0) in bottom left
if(MAT5)
  axis xy
end
