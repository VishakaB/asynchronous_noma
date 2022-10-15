SNRdb = 20;

SNR = 10 ^ (SNRdb/10);
n_t = ;      % number of Tx antennas
n_r = ;      % number of Rx antennas
rateth = ;             % rate threshold

%channel coefficients
H  = (randn(n_r, n_t) + sqrt(-1) * randn(n_r, n_t)) / sqrt(2);

%Calculate achievable rates


%Find average of achievable rates


%Check for outage



