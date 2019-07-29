function [T,F,P,O,A] = align_t_f(X,Y,Tmax)
% [T,F,P,O,A] = align_t_f(X,Y,Tmax)
%   X is a reference signal.  Y is a post-distortion signal which
%   has been (a) shifted in time by up to Tmax samples (b) shifted
%   in frequency (e.g. by demodulation on a single-side-band system
%   with a mistuned carrier).
%   The algorithm forms analytic versions of X and Y, multiplies
%   them at all possible relative skews out to Tmax points, then
%   takes the FFT of each product.  The global highest peak of this 
%   <tskew,FFTbin> matrix indicates the best time alignment and the 
%   best frequency shift to make Y match X.  O is returned as this
%   best-aligned version. T is the time skew, F is the frequency
%   shift, and P is the optimal phase of a sinusoid at F.
% 2011-04-22 Dan Ellis dpwe@ee.columbia.edu

% Make analytic versions of the signals
lx = length(X);
Xa = ifft(([1:lx]'<(lx/2)).*fft(X));
Ya = ifft(([1:lx]'<(lx/2)).*fft(Y));

%M = zeros(lx,2*Tmax+1);
%
%for ts = -Tmax:Tmax
%  M(:,ts+Tmax+1) = (fft(rotr(Ya,ts).*conj(Xa)));
%end
%% locate highest peak
%[vv,T] = max(max(abs(M)));
%[vv,F] = max(max(abs(M')));
%
%% When projecting back to the real axis, use the phase of this peak
%% component (i.e. the phase of the largest sinsoid we get from
%% conjugate multiplication, to rotate it back to the real axis)
%P = angle(M(F,T));
%T = T - Tmax - 1;

nfbins = lx;

if nargout > 4
  A = zeros(nfbins, 2*Tmax+1);
end


maxval = -Inf;
for ts = -Tmax:Tmax
  M = fft(rotr(Ya,ts).*conj(Xa));
  [vv,xx] = max(abs(M));
  if vv > maxval
    T = ts;
    F = xx;
    P = angle(M(xx));
    maxval = vv;
  end
  if nargout > 4
    A(:,ts+Tmax+1) = abs(M(1:nfbins));
  end
end


F = F - 1;
if F>(lx/2)
  F = F-lx;
end

% make the shifted version
O = real(rotr(Ya,T).*exp(-j*[0:lx-1]'/lx*F*2*pi-j*P));
  

  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = rotr(x, n)
% Rotate rows of x by n points
s = length(x);
while n < 0
  n = n + s;
end
y = [x((n+1):end); x(1:n)];
