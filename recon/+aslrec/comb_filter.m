function out = comb_filter(in,Nperiod, BW)
% function out = comb_filter(in,Nperiod, BW)
%
%  Making a comb filter as a bunch of notch filters concatenated
% it's meant to remove preiodic noise with harmonics
%
% inputs:
%   in = data vector
%   Nperiod = Nperiodicity of the repetition (in samples)
%   Q  = bandwith of the teeth in the comb function
%       
% outputs:
%   out = the cleaned up data
%

% make casorati matrix from 4D data (Nframes x Npix)
dims = size(in);
Npix = prod(dims(1:3));
Nframes = dims(4);
in = reshape(permute(in,[4 1:3]),[],Npix);


% calculate the frequency band of the main notch
% and the num. of harmonics to remove
w0 = Nperiod/(Nframes/2);
Nharmonics = floor(Nframes/2/Nperiod);

% zero padding
in = repmat(in,3,1);
% create and apply a notch filter for each harmonic
out= zeros(size(in));
parfor p=1:Npix
    out(:,p) = notchFilter(in(:,p), 1, w0, BW, Nharmonics  );
end
% remove padding
out = out(Nframes+1:end-Nframes, :);

% reshape data back to 4D matrix
out = reshape(out',dims);

end

%% old code to make a frq. response byt hand
%{
% create the frequency response of the fiter
% as a brute force comb, smoothed with a gaussian
f = linspace(0,1,Nframes);
% figure out the principal frequency
f0 = 2*Nframes/Nperiod
G = zeros(Nframes,1);
G(ceil(f0:f0:Nframes/2-2*Q)) = 1;

% gaussian smoothing of the comb determines the bandwidth of the teeth
s = exp(-linspace(-10,10, 20).^2 / Q^2);
G = conv(G,s);

% filter suppresses the desired bands, so flip it.
G = 1- G;
% clip the end - convolution byproduct.
G = G(ceil(length(s)/2):end-floor(length(s)/2));
% make is symmetric
G(end/2+1:end) = G(end/2:-1:1);
% show the response from 0 to 2*PI
plot(f,G);

% Calculate frequency spectrum of the input signals
IN = fft(in, [],1);
% apply the filter in frequency domain
OUT = IN .* G ;
% Calculate time domain of the output
out = ifft(OUT, [],1);
%}
