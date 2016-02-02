function [ x_wave ] = wavethresh_2D( x, waveletStages, waveletFilterName, wavecoeffS, lambdaWave, L )
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

    x_wave_helper = wavedec2(x,waveletStages,waveletFilterName);
    x_wave_helper = softthresh_real(x_wave_helper,lambdaWave/L);
    x_wave = waverec2(x_wave_helper,wavecoeffS,waveletFilterName);

end

