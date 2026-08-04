function [out] = comb_filter(in,period)

Fs = 1/length(in)            % Sampling frequency (relative)
f0 = Fs/period            % Fundamental interference frequency (relative)
Q = 35;             % Quality factor (higher = narrower notch)

x = in;

out = x;

wo = f0/(Fs/2)      % Normalized frequency
bw = wo/Q           % Bandwidth

% Filter each harmonic
for k = 1:floor((Fs)/f0 -1)
    fk = k*f0

    wo = fk/Fs      % Normalized frequency

    [b,a] = iirnotch(wo,bw);
    out = filtfilt(b,a,out); % Zero-phase filtering
    plot(abs(fft(out)))
    pause(0.2)
end

end