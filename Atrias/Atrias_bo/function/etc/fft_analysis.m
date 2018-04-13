function fft_analysis(t, signal, t_start, t_end, sample_rate)

    if nargin < 5
      sample_rate = 1000;
    end

    Fs = sample_rate;
    range = t>=t_start & t<=t_end;
    t_shortened = t(range);
    y = signal(range);
    L = length(t_shortened);


    figure('Name','FFT Analysis')
    clf

    subplot(211)

    plot(Fs*t_shortened,y)
    title('Signal over Time')
    xlabel('time (milliseconds)')


    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(y,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);

    % Plot single-sided amplitude spectrum.
    subplot(212)
    plot(f,2*abs(Y(1:NFFT/2+1))) 
    title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')

end