% Define simulation parameters
num_bits = 1e6; % Number of bits
snr_range = 0:2:30; % SNR range from 0 to 30 dB with 2 dB steps
signal_power = 1; % Assuming unit power for the signal

% Preallocate variables for speed
BER_manual = zeros(1, length(snr_range));

% Generate random binary data vector
bits = randi([0 1], 1, num_bits);

% Loop over the SNR range
for i = 1:length(snr_range)
    snr = snr_range(i);
    noise_power = signal_power / (10^(snr/10));
    noise = sqrt(noise_power) * randn(size(bits));
    Rx_sequence = bits + noise;

    % Decide whether the Rx_sequence is '1' or '0'
    detected_bits = Rx_sequence > 0.5;
    
    % Calculate the number of errors
    num_errors = sum(bits ~= detected_bits);
    
    % Save the probability of error for this SNR
    BER_manual(i) = num_errors / num_bits;
end

% Plot the BER curve against SNR for manual noise addition
figure;
semilogy(snr_range, BER_manual, 'b*-');
xlabel('Signal to Noise Ratio (SNR) [dB]');
ylabel('Bit Error Rate (BER)');
title('BER Curve against SNR (Manual Noise Addition)');
grid on;

 % Define simulation parameters (same as before)

% Preallocate variables for speed
BER_awgn = zeros(1, length(snr_range));

% Generate random binary data vector (same as before)

% Loop over the SNR range
for i = 1:length(snr_range)
    snr = snr_range(i);
    Rx_sequence = awgn(bits, snr, 'measured');

    % Decide whether the Rx_sequence is '1' or '0' (same as before)
        detected_bits = Rx_sequence > 0.5;
    % Calculate the number of errors (same as before)
        num_errors = sum(bits ~= detected_bits);
    % Save the probability of error for this SNR (same as before)
    BER_awgn(i) = num_errors / num_bits;
end

% Plot the BER curve against SNR for AWGN function
figure;
semilogy(snr_range, BER_awgn, 'r*-');
xlabel('Signal to Noise Ratio (SNR) [dB]');
ylabel('Bit Error Rate (BER)');
title('BER Curve against SNR (AWGN Function)');
grid on;

awgn requires Communications Toolbox.
 
