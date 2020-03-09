function [channel_statistics] = newRayleighChan(N, fm)
    
    % [channel_statistics] = newRayleighChan(N, fm)
    % [Usage]
    %   newRayleighChan creates channel coefficients that mimicks a
    %   flat-flading Rayleigh channel from inputs 'N' (number of sample
    %   points) and 'fm' (maximum Doppler shift frequency). 
    %   N should typically be a power of two.
    % >> procedure reference: Rappaport textbook Ch 5: pg 222
    
    N = N/2;             % compensate for IFFT for 2*N sample points
    % df = 2*fm / (N-1); % frequency spacing between adjacent spectral lines
    % T = 1/df;          % time duration of a fading waveform

    % column vectors
    randnoise1_pos = randn([N/2,1]) + 1j*randn([N/2,1]);
    randnoise1_neg = flipud(conj(randnoise1_pos));
    randnoise1 = [randnoise1_neg; randnoise1_pos];
    
    randnoise2_pos = randn([N/2,1]) + 1j*randn([N/2,1]);
    randdnoise2_neg = flipud(conj(randnoise2_pos));
    randnoise2 = [randdnoise2_neg; randnoise2_pos];

    f = linspace(-fm, fm, N);   % freq span vector
    fc = 0;                     % carrier frequency at 0 Hz
    SEz = 1.5 ./ ( pi*fm*sqrt( 1-( (f-fc)/fm ).^2 ) );
    
    % adjustment for 'Inf' at boundary of this Doppler Spectrum
    SEz(1) = 0;
    SEz(end) = 0;
    
    sqrtSEz = sqrt(SEz);
    randnoise1 = sqrtSEz .* randnoise1.';
    randnoise2 = sqrtSEz .* randnoise2.';
    
    time_randnoise1 = ifft(randnoise1, 2*N);
    time_randnoise2 = ifft(randnoise2, 2*N);
    
    % add the squares of each signal point in time to create an N-point time series
    sum = (time_randnoise1).^2 + (time_randnoise2).^2;
    channel_statistics = sqrt(sum).';
    
    % uncomment the following for Rayleigh (time) and Doppler Spectrum (freq) plot
    % figure;
    % n = linspace(0,length(channel_statistics),N*2);
    % plot(n(1:200)',20*log10(abs(channel_statistics(1:200))));
    % title('Rayleigh Fading Channel');
    % xlabel('Time Sample','FontWeight','bold','FontSize',12); 
    % ylabel('Amplitude','FontWeight','bold','FontSize',12);
    % 
    % figure;
    % plot(f,sqrtSEz);
    % title(sprintf('Doppler Spectrum for f_m=%d Hz, f_c = %d Hz',fm,fc));
    % xlabel('Frequency (Hz)','FontWeight','bold','FontSize',12);
    % ylabel('S_E_Z (f)','FontWeight','bold','FontSize',12);
end