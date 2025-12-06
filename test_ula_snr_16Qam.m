function test_ula_snr_v5
%% Parametri generali
%close all;

%% Direzioni
B = 6;
grad = -20; %20 correct, -25 error of -5 deg

theta_ref = deg2rad(0);
phi_ref   = deg2rad(-25); % destinatario

theta_int = deg2rad(0);
phi_int   = deg2rad(25);  % interferente

theta_int2 = deg2rad(0);
phi_int2   = deg2rad(-50);  % interferente2

%% Constellazione 16-QAM 
M_qam = 16;
qamSymbols = qammod(0:M_qam-1, M_qam, 'gray', 'UnitAveragePower', true);
k_qam = log2(M_qam);

%% ULA param
ulaSize = 3;
f0 = 1575.42e6;
c  = 299792458;
lambda = c/f0;
spacing = lambda/2;
antennaLocation = zeros(3,ulaSize);
antennaLocation(2,:) = (0:ulaSize-1) * spacing;

%% Maschera DM
azimuthWidth_deg   = 1;
amplitudeTolerance = 0.8;  % dB
phaseTolerance_deg = 20;

%% Colori simboli 16-QAM
myColor = lines(M_qam);

%% Simulazione pattern e scatter plot

N_pattern = 0; %1500;                      
tx_pattern = qamSymbols(randi(M_qam,N_pattern,1));  

rx_ref_pattern = complex(zeros(N_pattern,1));
rx_int_pattern = complex(zeros(N_pattern,1));

W_used = complex(zeros(ulaSize,N_pattern));
colors = zeros(N_pattern, 3);

% Prepara figura 
fh = figure(1);
fh.WindowState = 'maximized';
sgtitle('Directional Modulation con ULA e 16-QAM', ...
        'Interpreter','Latex','FontSize',20);

subplot(2,2,1); hold on; grid on; box on; axis([-90 +90 -30 0]);
xlabel('$\varphi$ [deg]','Interpreter','Latex');
ylabel('Ampiezza [dB]','Interpreter','Latex');

subplot(2,2,3); hold on; grid on; box on; axis([-90 +90 -180 +180]);
xlabel('$\varphi$ [deg]','Interpreter','Latex');
ylabel('Fase [deg]','Interpreter','Latex');

subplot(2,2,2); hold on; grid on; box on; axis square;
xlabel('I','Interpreter','Latex'); ylabel('Q','Interpreter','Latex');

subplot(2,2,4); hold on; grid on; box on; axis square;
xlabel('I','Interpreter','Latex'); ylabel('Q','Interpreter','Latex');

%% calcolo dei vettori
parfor n = 1:N_pattern
    [~,idx] = min(abs(tx_pattern(n) - qamSymbols));
    colors(n,:) = myColor(idx,:);
    
    %fprintf("Trovo vettore: %d/%d\n", n, N_pattern);

    [w, ok] = FindSteeringVector(ulaSize, lambda, antennaLocation, ...
        theta_ref, phi_ref, azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg);

    if ~ok
        warning('Nessuna soluzione trovata al simbolo %d\n', n);
        continue;
    end

    W_used(:,n) = w;

    a_ref = exp(1i*myWaveVector(lambda,theta_ref,phi_ref)'*antennaLocation).';
    a_int = exp(1i*myWaveVector(lambda,theta_int,phi_int)'*antennaLocation).';

    rx_ref_pattern(n) = (w'*a_ref)/ulaSize * tx_pattern(n);
    rx_int_pattern(n) = (w'*a_int)/ulaSize * tx_pattern(n);
end

%% PLOT scatter plot 
subplot(2,2,2); hold on;
for n = 1:N_pattern
    plot(real(rx_ref_pattern(n)), imag(rx_ref_pattern(n)), 'o', ...
         'MarkerSize',6, 'MarkerFaceColor', colors(n,:), 'MarkerEdgeColor','k');
end
%title('Ricevuto - destinatario (I/Q)');

subplot(2,2,4); hold on;
for n = 1:N_pattern
    plot(real(rx_int_pattern(n)), imag(rx_int_pattern(n)), 'o', ...
         'MarkerSize',6, 'MarkerFaceColor', colors(n,:), 'MarkerEdgeColor','k');
end
%title('Ricevuto - interferente (I/Q)');

%% PLOT ampiezza e fase array per ogni simbolo
phi = linspace(-pi/2, pi/2, 1001);
theta = 0;

for n = 1:N_pattern
    [~, idxSym] = min(abs(tx_pattern(n) - qamSymbols));
    thisColor   = myColor(idxSym,:);
    w = W_used(:,n);
    if all(w==0), continue; end

    arrayResponse = zeros(1,length(phi));
    for idx = 1:length(phi)
        arrayResponse(idx) = (w'*exp(1i*myWaveVector(lambda,theta,phi(idx))' ...
            *antennaLocation).') / ulaSize;
    end

    % Ampiezza dB
    subplot(2,2,1);
    plot(rad2deg(phi), pow2db(abs(arrayResponse)), '-', ...
         'LineWidth', 0.8, 'Color', thisColor);

    % Fase 
    subplot(2,2,3);
    plot(rad2deg(phi), rad2deg(angle(arrayResponse)), '-', ...
         'LineWidth', 0.8, 'Color', thisColor);
end

%% === PULSANTE PER SALVATAGGIO ===
%uicontrol('Style', 'pushbutton', 'String', 'Salva Vettori',...
%    'Units','normalized','Position', [0.8 0.92 0.15 0.05], ...
%    'FontSize',12, 'Callback', @(src,event) salvaVettori(W_used,ulaSize));

%% ===============================================================
%   Simulazione BER 16-QAM al variare di Eb/N0 usando i vettori salvati
% ===============================================================

fileName = sprintf('steeringVector%d_%d-%d_%d_%d_%d.mat', ulaSize, grad, azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg, B);
if isfile(fileName)
    load(fileName, 'allVectors');
    nVectors = numel(allVectors);
    fprintf("Sono stati trovati %d set di vettori salvati per ULA=%d.\n", nVectors, ulaSize);
else
    error("File %s non trovato!", fileName);
end

EbN0_dB     = 0:1:15;       
nMax_error  = 1000;          
BER_dest_sim = zeros(size(EbN0_dB));
BER_int_sim  = zeros(size(EbN0_dB));
BER_int2_sim = zeros(size(EbN0_dB));

%% Loop Eb/N0
for k = 1:length(EbN0_dB)
    SNR_lin = 10^(EbN0_dB(k)/10);
    
    % Inizializzazione contatori
    N_symbols     = 0; 
    N_errors_dest = 0; 
    N_errors_int  = 0; 
    N_errors_int2 = 0;

    % Simulazione fino a nMax_error
    while N_errors_dest < nMax_error

        % Seleziona un vettore casuale dal file
        idxVec = randi(nVectors);
        W_test = allVectors{idxVec};
        col = randi(size(W_test,2));
        w = W_test(:,col);
        %caso ottimo ula = 5
        %w = [1.0000 + 0.0000i; 0.4762 - 0.8793i; -0.5465 - 0.8375i; -0.9967 + 0.0818i; -0.4027 + 0.9153i];

        %steering vectors
        a_ref  = exp(1i*myWaveVector(lambda,theta_ref,phi_ref)'*antennaLocation).';
        a_int  = exp(1i*myWaveVector(lambda,theta_int,phi_int)'*antennaLocation).';
        a_int2 = exp(1i*myWaveVector(lambda,theta_int2,phi_int2)'*antennaLocation).';

        % simbolo 16QAM casuale
        symIdx = randi(M_qam);
        txSym  = qamSymbols(symIdx);
        qamBits = de2bi(0:M_qam-1, log2(M_qam), 'left-msb');
        txBits = qamBits(symIdx,:);

        gain_dest = abs((w'*a_ref)/ulaSize);
        gain_int  = abs((w'*a_int)/ulaSize);
        gain_int2 = abs((w'*a_int2)/ulaSize);
        
        % normalizzazione 
        rx_dest = gain_dest  * txSym;  
        rx_int  = gain_int   * txSym;
        rx_int2 = gain_int2  * txSym;
        
        % rumore
        sigma2 = 1 / (SNR_lin * k_qam);
        
        noise_dest = sqrt(sigma2 / 2) * (randn + 1i*randn);
        noise_int  = sqrt(sigma2 / 2) * (randn + 1i*randn);
        noise_int2 = sqrt(sigma2 / 2) * (randn + 1i*randn);
        
        % aggiunta rumore al segnale
        rx_dest_rumore = rx_dest + noise_dest;
        rx_int_rumore  = rx_int  + noise_int;
        rx_int2_rumore = rx_int2 + noise_int2;

        % decisione ML
        [~, idx] = min(abs(rx_dest_rumore - qamSymbols));
        rxBits_dest = qamBits(idx,:);

        [~, idx] = min(abs(rx_int_rumore - qamSymbols));
        rxBits_int = qamBits(idx,:);

        [~, idx] = min(abs(rx_int2_rumore - qamSymbols));
        rxBits_int2 = qamBits(idx,:);

        % aggiorno contatori
        N_symbols = N_symbols + 1;
        N_errors_dest = N_errors_dest + sum(txBits ~= rxBits_dest);
        N_errors_int  = N_errors_int  + sum(txBits ~= rxBits_int);
        N_errors_int2 = N_errors_int2 + sum(txBits ~= rxBits_int2);
    end

    BER_dest_sim(k) = N_errors_dest / (k_qam * N_symbols);
    BER_int_sim(k)  = N_errors_int  / (k_qam * N_symbols);
    BER_int2_sim(k) = N_errors_int2 / (k_qam * N_symbols);

    fprintf("Eb/N0 = %d dB | Simboli: %d | BER_dest = %.4e\n", ...
        EbN0_dB(k), N_symbols, BER_dest_sim(k));
end

%% BER teorica 16-QAM
EbN0_lin = 10.^(EbN0_dB/10);
BER_theory  = (1-(1-2*(1-1/sqrt(M_qam))*qfunc(sqrt(3*k_qam*EbN0_lin/(M_qam-1)))).^2)/k_qam;

%% Plot BER 
figure;
semilogy(EbN0_dB, BER_theory, 'b--','LineWidth',2); hold on;
semilogy(EbN0_dB, BER_dest_sim, 'g-','LineWidth',2);
semilogy(EbN0_dB, BER_int_sim, 'r--o','LineWidth',2, 'MarkerSize',6);
semilogy(EbN0_dB, BER_int2_sim, 'Color', [1 0.5 0], 'LineStyle','-.', 'LineWidth',2, 'Marker','^', 'MarkerSize',6);

grid on;
ax = gca;
ax.FontSize=12;

xlabel('$E_b\,/N_0$ [dB]', 'Interpreter', 'latex', 'FontSize',20);
ylabel('BER', 'Interpreter', 'latex', 'FontSize',20);

%title(sprintf('BER QAM, ULA\\_SIZE=%d', ulaSize), 'Interpreter', 'tex');
legend('Theoretical','Intended receiver','Eavesdropper 1','Eavesdropper 2', 'Location','southwest', 'Interpreter','latex', 'FontSize', 16);

end

%% === Funzioni di supporto ===
function k = myWaveVector(lambda, theta, phi)
    k = 2*pi/lambda * [cos(theta).*cos(phi); cos(theta).*sin(phi); sin(theta)];
end

function [steeringVector, solutionIsFound] = FindSteeringVector(ulaSize, lambda, antennaLocation, ...
    theta_ref, phi_ref, azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg)

    maxIter = 1e7; 
    phi = linspace(-pi/2,pi/2,201);
    theta = 0;
    solutionIsFound = false;
    steeringVector = complex(zeros(ulaSize,1));

    for n = 1:maxIter
        B = 6;
        numLevels = 2^B;
        levels = (0:numLevels-1)'*(2*pi/numLevels);
        idx = randi(numLevels, ulaSize,1);
        phaseVector = levels(idx) + deg2rad(2)*randn(ulaSize,1);
        steeringVector = exp(1i*phaseVector);

        arrayResponse = zeros(1,length(phi));
        for index = 1:length(phi)
            arrayResponse(index) = (steeringVector' * ...
                exp(1i*myWaveVector(lambda,theta,phi(index))'*antennaLocation).') / ulaSize;
        end

        sel = (phi >= phi_ref-deg2rad(azimuthWidth_deg)/2) & ...
              (phi <= phi_ref+deg2rad(azimuthWidth_deg)/2);
        if ~any(sel), continue; end
        selectedArrayResponse = arrayResponse(sel);

        amplitudes = abs(selectedArrayResponse);
        phases_deg = rad2deg(angle(selectedArrayResponse));
        maxAmp_dB = 20*log10(max(amplitudes));
        minAmp_dB = 20*log10(min(amplitudes));

        if maxAmp_dB <= 0 && minAmp_dB >= -amplitudeTolerance && ...
           max(phases_deg) <= phaseTolerance_deg/2 && ...
           min(phases_deg) >= -phaseTolerance_deg/2
            solutionIsFound = true;
            return;
        end
    end
end

function salvaVettori(W_used, ulaSize)
    fileName = sprintf('steeringVector%d.mat', ulaSize);
    if isfile(fileName)
        S = load(fileName, 'allVectors');
        if isfield(S,'allVectors')
            allVectors = S.allVectors;
            if ~iscell(allVectors)
                allVectors = {allVectors};
            end
        else
            allVectors = {};
        end
    else
        allVectors = {};
    end
    newVectors = num2cell(W_used, 1);
    newVectors = newVectors(:);
    allVectors = [allVectors; newVectors];
    save(fileName, 'allVectors', '-v7.3');
    fprintf('Vettori salvati in %s (totale set: %d)\n', fileName, numel(allVectors));
end
