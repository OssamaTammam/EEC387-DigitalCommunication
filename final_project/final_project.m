% Simulation Parameters
N = 1e5;  % Number of bits
SNR_dB = 0:4:60;  % SNR range in dB
SNR = 10.^(SNR_dB/10);  % SNR in linear scale

% Generate random binary data vector
data = randi([0 1], N, 1);

% Initialize BER arrays
BER_OOK = zeros(length(SNR_dB), 1);
BER_PRK = zeros(length(SNR_dB), 1);
BER_FSK = zeros(length(SNR_dB), 1);

% Loop over different SNR values
for i = 1:length(SNR_dB)
    % OOK Modulation and Demodulation
    OOK_signal = data;  % OOK modulation: no change
    noise_OOK = sqrt(1/(2*SNR(i)))*randn(N, 1);
    OOK_received = OOK_signal + noise_OOK;  % Add noise
    OOK_demod = OOK_received > 0.5;  % Demodulation
    BER_OOK(i) = sum(data ~= OOK_demod) / N;  % Calculate BER

    % PRK (BPSK) Modulation and Demodulation
    PRK_signal = 2*data - 1;  % PRK modulation: 0 -> -1, 1 -> 1
    noise_PRK = sqrt(1/(2*SNR(i)))*randn(N, 1);
    PRK_received = PRK_signal + noise_PRK;  % Add noise
    PRK_demod = PRK_received > 0;  % Demodulation
    BER_PRK(i) = sum(data ~= PRK_demod) / N;  % Calculate BER

    % Orthogonal-FSK Modulation and Demodulation
    FSK_signal = (data == 0) + 1i*(data == 1);  % FSK modulation
    noise_FSK = sqrt(1/(2*SNR(i)))*(randn(N, 1) + 1i*randn(N, 1));
    FSK_received = FSK_signal + noise_FSK;  % Add noise
    FSK_demod = real(FSK_received) < imag(FSK_received);  % Demodulation
    BER_FSK(i) = sum(data ~= FSK_demod) / N;  % Calculate BER
end

% Plot results
figure;
semilogy(SNR_dB, BER_OOK, '-o', 'DisplayName', 'OOK');
hold on;
semilogy(SNR_dB, BER_PRK, '-x', 'DisplayName', 'PRK');
semilogy(SNR_dB, BER_FSK, '-s', 'DisplayName', 'Orthogonal-FSK');
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR for Different Modulation Schemes');
legend;
grid on;

% Evaluate using MATLAB built-in functions
mod_types = {'pam', 'psk', 'fsk'};
BER_builtin = zeros(length(SNR_dB), length(mod_types));

for j = 1:length(mod_types)
    for i = 1:length(SNR_dB)
        M = 2;  % Binary modulation
        switch mod_types{j}
            case 'pam'
                hMod = comm.PAMModulator(M);
                hDemod = comm.PAMDemodulator(M);
            case 'psk'
                hMod = comm.PSKModulator(M, 'BitInput', true);
                hDemod = comm.PSKDemodulator(M, 'BitOutput', true);
            case 'fsk'
                hMod = comm.FSKModulator(M, 'BitInput', true);
                hDemod = comm.FSKDemodulator(M, 'BitOutput', true);
        end
        y = step(hMod, data);
        y_noisy = awgn(y, SNR_dB(i), 'measured');
        z = step(hDemod, y_noisy);
        BER_builtin(i, j) = biterr(data, z) / N;
    end
end

% Plot built-in function results
figure;
semilogy(SNR_dB, BER_builtin(:, 1), '-o', 'DisplayName', 'Built-in PAM');
hold on;
semilogy(SNR_dB, BER_builtin(:, 2), '-x', 'DisplayName', 'Built-in PSK');
semilogy(SNR_dB, BER_builtin(:, 3), '-s', 'DisplayName', 'Built-in FSK');
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR for Different Modulation Schemes (Built-in)');
legend;
grid on;

% 16QAM Modulation using qammod and qamdemod
M = 16;
dataSymbols = randi([0 M-1], N, 1);  % Generate random symbols for 16QAM
BER_16QAM = zeros(length(SNR_dB), 1);

for i = 1:length(SNR_dB)
    % Modulate
    y = qammod(dataSymbols, M, 'UnitAveragePower', true);
    % Add noise
    y_noisy = awgn(y, SNR_dB(i), 'measured');
    % Demodulate
    z = qamdemod(y_noisy, M, 'UnitAveragePower', true);
    % Calculate BER
    [~, BER_16QAM(i)] = biterr(dataSymbols, z);
end

% Plot 16QAM results
figure;
semilogy(SNR_dB, BER_16QAM, '-o', 'DisplayName', '16QAM');
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR for 16QAM');
legend;
grid on;
