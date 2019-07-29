function S = nist_stnr_m(D,SR,doplot)
% S = nist_stnr(D,SR,doplot)
%    Calculate NIST STNR actually in Matlab, attempting to duplicate
%    stnr -c algorithm
% 2011-08-02 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; doplot = 0; end

verbose = 0;

% constants from snr.h
MILI_SEC = 20.0;
PEAK_LEVEL = 0.95;

BINS = 500;
%SMOOTH_BINS = 7;
SMOOTH_BINS = 15;
CODEC_SMOOTH_BINS = 15;
LOW = -28.125;
HIGH = 96.875;

BLOCKSIZE = 2048;

% from hist.h
PEAK = 0;
TROUGH = 1;
PEAK_WIDTH = 3;

% algorithm is to form a power histogram over 20ms windows

frame_width = SR / 1000 * MILI_SEC;
frame_adv = frame_width / 2;

% calculate power, assuming short samples
D2 = (D*16384).^2;

nhops = floor(length(D2)/frame_adv);
D2 = reshape(D2(1:(nhops*frame_adv)), frame_adv, nhops);
% Power in each half window
P2 = mean(D2);
% Power in overlapping windows
Pdb = 10*log10(conv([1 1],P2));

% Histogram
hvals = LOW:(HIGH-LOW)/BINS:HIGH;
power_hist = hist(Pdb, hvals);

% stnr -c algorithm
unspiked_hist = medianf(power_hist,3);
presmooth_hist = conv2(unspiked_hist, ...
                       ones(1,CODEC_SMOOTH_BINS*2+1)/(CODEC_SMOOTH_BINS*2+1), ...
                       'same');
smoothed_hist = conv2(presmooth_hist, ...
                      ones(1,SMOOTH_BINS*2+1)/(SMOOTH_BINS*2+1), ...
                      'same');

if doplot
  subplot(411)
  plot(hvals,power_hist);
  title('power_hist');
  subplot(412)
  plot(hvals,unspiked_hist);
  title('unspiked_hist');
  subplot(413)
  plot(hvals,presmooth_hist);
  title('presmooth_hist');
  subplot(414)
  plot(hvals,smoothed_hist);
  title('smoothed_hist');
end

% assume to begin with that we don't find any extrema */
first_peak=BINS;
first_trough=BINS;
second_peak=BINS;
second_trough=BINS;

max_val = max(smoothed_hist);

% now look for the extrema, sequentially */

% find the noise peak; it should be reasonably big */
starting_point=0;
first_peak = locate_extremum(smoothed_hist,starting_point,BINS,PEAK);
while (10*smoothed_hist(1+first_peak)) < max_val
    first_peak = locate_extremum(smoothed_hist,starting_point,BINS,PEAK);
    starting_point = first_peak+1;
end

% now find the rest */
first_trough = locate_extremum(smoothed_hist,first_peak+1,BINS,TROUGH);
second_peak = locate_extremum(smoothed_hist,first_trough+1,BINS,PEAK);
second_trough = locate_extremum(smoothed_hist,second_peak+1,BINS,TROUGH);

if verbose
  fprintf(1, ...
          'peak=%d (%5.2f) trough=%d (%5.2f) peak=%d (%5.2f) trough=%d (%5.2f)\n', ...
          first_peak, ...
          pick_center(smoothed_hist,first_peak), ...
          first_trough, ...
          pick_center(smoothed_hist,first_trough), ...
          second_peak, ...
          pick_center(smoothed_hist,second_peak), ...
          second_trough, ...
          pick_center(smoothed_hist,second_trough));
end

if first_peak==BINS
  if verbose
    fprintf(1, ['I can''t find the first peak of the power distribution. ' ...
                'Is this a null file?\n']);
  end
  S = 0;
  return
end

noise_lvl = pick_center(smoothed_hist, first_peak);

if first_trough==BINS
  if verbose
    fprintf(1, ['Can''t find first trough. I''ll do my best from ' ...
                'here...\n']);
  end
  
  for i=0:(first_peak-1)
    full_hist(1+i) = 0;
  end
  cross_lvl = -Inf;
  speech_lvl = percentile_hist(unspiked_hist,BINS,PEAK_LEVEL);
  S = speech_lvl - noise_lvl;
  return;
end

for i=0:first_trough-1
  unspiked_hist(1+i)=0;
end

if second_peak==BINS
  if verbose
    fprintf(1,'Can''t find second peak.');
  end
  
  cross_lvl = -Inf;
  speech_lvl = percentile_hist(unspiked_hist,BINS,PEAK_LEVEL);
  S = speech_lvl - noise_lvl;
  return;
end

% check for bogus hump */
if 60*(smoothed_hist(1+second_peak)-smoothed_hist(1+first_trough)) ...
      < smoothed_hist(1+first_peak)
  cross_lvl = -Inf;
else
  cross_lvl = pick_center(smoothed_hist, second_peak);
end

if second_trough == BINS
  second_lim = second_peak;
else
  second_lim = second_trough;
end

for i = 0:second_lim-1
  unspiked_hist(1+i)=0;
end

speech_lvl = percentile_hist(unspiked_hist,BINS,PEAK_LEVEL);

S = speech_lvl - noise_lvl;

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = locate_extremum(h,from,to,type)
% from hist.c

% from hist.h
PEAK = 0;
TROUGH = 1;
PEAK_WIDTH = 3;

for i= from+PEAK_WIDTH:to-PEAK_WIDTH-1
  if h(1+i) == 0  % not interested in extrema at 0 
    continue;
  end
  
  extremum=1; % assume it's an extremum to begin with
  pre_swing=0;
  post_swing=0;
  swing_loc=i-PEAK_WIDTH;
  for j = i-PEAK_WIDTH:i-1 % check the preceding samples
    if type==PEAK
      if h(1+j) > h(1+j+1)
        extremum=0;
        break;
      end
      if h(1+j) ~= h(1+j+1)
        pre_swing=1;
      end
    else % type == TROUGH */
      if h(1+j) < h(1+j+1)
        extremum=0;
        break;
      end
      if h(1+j) ~= h(1+j+1)
        pre_swing=1;
      end
    end
  end
      
  if extremum == 0
    continue; 
  end
    
  for j=i:i+PEAK_WIDTH-1  % check the subsequent samples
    if type==PEAK
      if h(1+j) < h(1+j+1)
        extremum=0;
        break;
      end
      if h(1+j) ~= h(1+j+1)
        post_swing=1;
      end
    else % type == TROUGH
      if h(1+j) > h(1+j+1)
        extremum=0;
        break;
      end
      if h(1+j) ~= h(1+j+1)
        post_swing=1;
      end
    end
  end

  % check to make sure it isn't a step function
  % this kind of check is necessary if the peak is wider than the window
  if (((pre_swing+post_swing)<=1)&&(extremum))
    for k=i:-1:(from+1)
      diff1 = h(1+k-1) - h(1+k);
      if diff1 ~= 0
        break;
      end
    end
    swing_loc=k;
    for k=i:to-1-1   % find next swing
      diff2 = h(1+k) - h(1+k+1);
      if diff2 ~=0
        break;
      end
    end
    next_swing_loc=k;
    if ((type==PEAK)&&((diff1>0)||(diff2<0))) continue; end  % no dice
    if ((type==TROUGH)&&((diff1<0)||(diff2>0))) continue; end % ditto 
    
    % otherwise, the peak is at the mid-point of this plateau
    retval = round((swing_loc+next_swing_loc)/2);
    return;
  end
  if (extremum)
    retval = i;
    return
  end
end

retval = to;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val =  pick_center(h, bin)
% from snr.c

BINS = 500;
LOW = -28.125;
HIGH = 96.875;

step = (HIGH-LOW)/BINS;

val = LOW+step*(bin+0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = percentile_hist(h, bins, percent)

cumm = cumsum(h);
bin = min(find(cumm >= percent*cumm(end)))-1;

val = pick_center(h, bin);
