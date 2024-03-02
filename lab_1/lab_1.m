% Clear the command window and workspace
clc;
clear all;

% Define parameters
N_bits = 10;
N_s_pb = 100;
bit_seq = randi([0, 1], 1, N_bits);
Rb = 1; % 1 bit per second
Tb = 1; % 1/Rb
Ts = Tb / N_bits; % Time for one sample

% Polar NRZ encoding
Polar_NRZ = zeros(1, N_bits * N_s_pb);
for i = 0:1:N_bits-1
    if (bit_seq(i+1) == 1)
        Polar_NRZ(i*N_s_pb+1:(i+1)*N_s_pb) = 1;
    else
        if (bit_seq(i+1) == 0)
            Polar_NRZ(i*N_s_pb+1:(i+1)*N_s_pb) = -1;
        end
    end
end

% Plot Polar NRZ
figure;
subplot(6, 1, 1);
plot(Polar_NRZ, 'LineWidth', 2); grid on;
axis([1 N_bits*Tb*N_s_pb -1.5 1.5]);
xlabel('Time (mSec)');
ylabel('Amplitude (Volts)');
title('Polar NRZ');

% Inverted Polar NRZ encoding
Polar_NRZ_I = zeros(1, N_bits * N_s_pb);
for i = 0:1:N_bits-1
    if (i == 0)
        Polar_NRZ_I(i*N_s_pb+1:(i+1)*N_s_pb) = 1;
    else
        if (bit_seq(i+1) == 1)
            if (Polar_NRZ_I(i*N_s_pb-1) == 1)
                Polar_NRZ_I(i*N_s_pb+1:(i+1)*N_s_pb) = -1;
            else
                Polar_NRZ_I(i*N_s_pb+1:(i+1)*N_s_pb) = 1;
            end
        else
            if (bit_seq(i+1) == 0)
                Polar_NRZ_I(i*N_s_pb+1:(i+1)*N_s_pb) = Polar_NRZ_I(i*N_s_pb-1);
            end
        end
    end
end

% Plot Inverted Polar NRZ
subplot(6, 1, 2);
plot(Polar_NRZ_I, 'LineWidth', 2); grid on;
axis([1 N_bits*Tb*N_s_pb -1.5 1.5]);
xlabel('Time (mSec)');
ylabel('Amplitude (Volts)');
title('Inverted Polar NRZ');

% Polar RZ encoding
N_s_pb_RZ = 100;
Polar_RZ = zeros(1, N_bits * N_s_pb_RZ * 2);
for i = 0:1:N_bits-1
    if (bit_seq(i+1) == 1)
        Polar_RZ(i*N_s_pb_RZ+1:i*N_s_pb_RZ+50) = 1;
        Polar_RZ(i*N_s_pb_RZ+1+50:(i+1)*N_s_pb_RZ) = 0;
    else
        if (bit_seq(i+1) == 0)
            Polar_RZ(i*N_s_pb_RZ+1:i*N_s_pb_RZ+50) = -1;
            Polar_RZ(i*N_s_pb_RZ+1+50:(i+1)*N_s_pb_RZ) = 0;
        end
    end
end

% Plot Polar RZ
subplot(6, 1, 3);
plot(Polar_RZ, 'LineWidth', 2); grid on;
axis([1 N_bits*Tb*N_s_pb -1.5 1.5]);
xlabel('Time (mSec)');
ylabel('Amplitude (Volts)');
title('Polar RZ');

% Bipolar (AMI) encoding
N_s_pb_BPZ = 100;
Cnt = 0;
BPZ = zeros(1, N_bits * N_s_pb_RZ * 2);
for i = 0:1:N_bits-1
    if (bit_seq(i+1) == 1)
        if (mod(Cnt, 2) == 0)
            BPZ(i*N_s_pb_BPZ+1:(i+1)*N_s_pb_BPZ) = 1;
            Cnt = Cnt + 1;
        else
            BPZ(i*N_s_pb_BPZ+1:(i+1)*N_s_pb_BPZ) = -1;
            Cnt = Cnt + 1;
        end
    else
        if (bit_seq(i+1) == 0)
            BPZ(i*N_s_pb_BPZ+1:(i+1)*N_s_pb_BPZ) = 0;
        end
    end
end

% Plot Bipolar (AMI)
subplot(6, 1, 4);
plot(BPZ, 'LineWidth', 2); grid on;
axis([1 N_bits*Tb*N_s_pb -1.5 1.5]);
xlabel('Time (mSec)');
ylabel('Amplitude (Volts)');
title('Bipolar (AMI)');

% Manchester encoding
N_s_pb_Mufc = 100;
MUFC = zeros(1, N_bits * N_s_pb_RZ * 2);
for i = 0:1:N_bits-1
    if (bit_seq(i+1) == 1)
        MUFC(i*N_s_pb_RZ+1:i*N_s_pb_RZ+50) = 1;
        MUFC(i*N_s_pb_RZ+1+50:(i+1)*N_s_pb_RZ) = -1;
    else
        if (bit_seq(i+1) == 0)
            MUFC(i*N_s_pb_RZ+1:i*N_s_pb_RZ+50) = -1;
            MUFC(i*N_s_pb_RZ+1+50:(i+1)*N_s_pb_RZ) = 1;
        end
    end
end

% Plot Manchester
subplot(6, 1, 5);
plot(MUFC, 'LineWidth', 2); grid on;
axis([1 N_bits*Tb*N_s_pb -1.5 1.5]);
xlabel('Time (mSec)');
ylabel('Amplitude (Volts)');
title('Manchester');

% MULT3 encoding
N_s_pb_MULT3 = 100;
States_cnt = 0;
MULT3 = zeros(1, N_bits * N_s_pb_RZ * 2);
for i = 0:1:N_bits-1
    if (bit_seq(i+1) == 1)
        if (States_cnt == 0)
            MULT3(i*N_s_pb_MULT3+1:(i+1)*N_s_pb_MULT3) = 1;
            States_cnt = States_cnt + 1;
        else
            if (States_cnt == 1)
                MULT3(i*N_s_pb_MULT3+1:(i+1)*N_s_pb_MULT3) = 0;
                States_cnt = States_cnt + 1;
            else
                if (States_cnt == 2)
                    MULT3(i*N_s_pb_MULT3+1:(i+1)*N_s_pb_MULT3) = -1;
                    States_cnt = States_cnt + 1;
                else
                    if (States_cnt == 3)
                        MULT3(i*N_s_pb_MULT3+1:(i+1)*N_s_pb_MULT3) = 0;
                        States_cnt = States_cnt + 1;
                    else
                        if (States_cnt == 4)
                            MULT3(i*N_s_pb_MULT3+1:(i+1)*N_s_pb_MULT3) = 1;
                            States_cnt = 0;
                        end
                    end
                end
            end
        end
    else
        if (bit_seq(i+1) == 0)
            if (i == 0)
                MULT3(i*N_s_pb_MULT3+1:(i+1)*N_s_pb_MULT3) = 0;
            else
                MULT3(i*N_s_pb_MULT3+1:(i+1)*N_s_pb_MULT3) = MULT3(i*N_s_pb_MULT3);
            end
        end
    end
end

% Plot MULT3
subplot(6, 1, 6);
plot(MULT3, 'LineWidth', 2); grid on;
axis([1 N_bits*Tb*N_s_pb -1.5 1.5]);
xlabel('Time (mSec)');
ylabel('Amplitude (Volts)');
title('MULT3');

% Plot power spectral densities
figure;
subplot(6, 1, 1);
T_B = 0.1;
f = linspace(0, 100, 10000);
PSD_POLAR_NRZ = T_B .* sinc(f .* T_B) .* sinc(f .* T_B);
xlim([0 100]);
plot(f, PSD_POLAR_NRZ);
xlabel('f');
ylabel('S(f)');
title('PSD POLAR NRZ');

subplot(6, 1, 2);
T_B = 0.1;
f = linspace(0, 100, 10000);
PSD_POLAR_NRZ_I = (T_B/4) .* sinc(f .* T_B) .* sinc(f .* T_B);
xlim([0 100]);
plot(f, PSD_POLAR_NRZ_I);
xlabel('f');
ylabel('S(f)');
title('PSD POLAR NRZ_I');

subplot(6, 1, 3);
T_B = 0.1;
T_B_2 = 0.05;
T_B_4 = 0.025;
f = linspace(0, 100, 10000);
PSD_POLAR_RZ = T_B_4 .* sinc(f .* T_B_2) .* sinc(f .* T_B_2);
plot(f, PSD_POLAR_RZ);
xlim([0 100]);
xlabel('f');
ylabel('S(f)');
title('PSD Polar RZ');

subplot(6, 1, 4);
T_B = 0.1;
T_B_2 = 0.05;
T_B_4 = 0.025;
f = linspace(0, 100, 10000);
PSD_BPOLAR_RZ = T_B_4 .* sinc(f .* T_B_2) .* sinc(f .* T_B_2) .* sin(pi .* f .* T_B) .* sin(pi .* f .* T_B);
plot(f, PSD_POLAR_RZ);
xlabel('f');
ylabel('S(f)');
title('PSD Bipolar RZ');

subplot(6, 1, 5);
T_B = 0.1;
T_B_2 = 0.05;
T_B_4 = 0.025;
f = linspace(0, 100, 10000);
PSD_MUFC = T_B .* sinc(f .* T_B_2) .* sinc(f .* T_B_2) .* sin(pi .* f .* T_B_2) .* sin(pi .* f .* T_B_2);
xlim([0 100]);
plot(f, PSD_MUFC);
xlabel('f');
ylabel('S(f)');
title('PSD MUFC');

subplot(6, 1, 6);
T_B = 0.1;
f = linspace(0, 100, 10000);
PSD_MLT3 = (4 / (T_B * T_B)) .* sin(pi .* f .* T_B/2) .* sin(pi .* f .* T_B/2) ./ ((pi .* f .* T_B) .* ((pi .* f .* T_B)));
xlim([0 100]);
plot(f, PSD_MLT3);
xlabel('f');
ylabel('S(f)');
title('PSD MLT3');
