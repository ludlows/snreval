# snreval

## Objective measures of speech quality/SNR

This collection of Matlab functions calculates a set of objective speech quality measures, mostly focused around some version of SNR (i.e. speech energy to nonspeech energy ratio). The measures are:

    **NIST STNR** - see http://labrosa.ee.columbia.edu/~dpwe/tmp/nist/doc/stnr.txt 

    **WADA SNR**  - see http://www.cs.cmu.edu/~robust/Papers/KimSternIS08.pdf

    **BSS_EVAL**  - see http://bass-db.gforge.inria.fr/bss_eval/ 

    **PESQ**      - see https://en.wikipedia.org/wiki/PESQ  

    **SNR_VAD**   - the "extra" energy in regions designated as speech by some kind of voice activity detection (VAD) when compared to the energy of the "gaps" in-between.


The code is based on the version on the link ( https://labrosa.ee.columbia.edu/projects/snreval ).

## install PESQ

see https://github.com/ludlows/pesq-mex , you need to run compile_test.m to make sure the compiling work if you didn't compile before.

PESQ_MEX_PATH = '~/Desktop/pesq-mex/'; % configure your pesq-mex path
addpath(PESQ_MEX_PATH)
addpath([PESQ_MEX_PATH, 'bin/'])


## Usage

    snreval('sample/arabic_400mhz.wav','-vad','sample/arabic_400mhz-vad.txt','-clean','sample/arabic_source.wav');

#==========================================
Target File: sample/arabic_400mhz.wav
 time range: 0.00000-60.00000 s
   Ref File: sample/arabic_source.wav
 Targ delay: -0.237 s
NIST STNR = 23.75000 dB
WADA SNR  = 19.62242 dB
   SNRvad = 12.34779 dB
      SAR = 9.07845 dB
 PESQ NARROW MOS-LQO = 2.27732
 PESQ NARROW MOS     = 2.61291
#==========================================


    snreval('sample/arabic_400mhz.wav','-guessvad',1); 

#==========================================
Target File: sample/arabic_400mhz.wav
 time range: 0.00000-60.00000 s
NIST STNR = 23.75000 dB
WADA SNR  = 19.62242 dB
   SNRvad = 8.10062 dB
#==========================================


    snreval('sample/arabic_400mhz.wav','-vad','sample/arabic_400mhz-vad.txt');

#==========================================
Target File: sample/arabic_400mhz.wav
 time range: 0.00000-60.00000 s
NIST STNR = 23.75000 dB
WADA SNR  = 19.62242 dB
   SNRvad = 12.34779 dB
#==========================================










