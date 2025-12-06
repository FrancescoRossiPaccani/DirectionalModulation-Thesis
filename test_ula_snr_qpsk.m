function test_ula_snr_QPSK
%% Parametri generali
%close all;

%% Direzioni

grad = -20;
B = 6;

theta_ref = deg2rad(0);
phi_ref   = deg2rad(-20); % destinatario

theta_int = deg2rad(0);
phi_int   = deg2rad(25);  % eavesdropper1

theta_int2 = deg2rad(0);
phi_int2   = deg2rad(-50);  % eavesdropper2

%% Constellazione QPSK
qpskConstellation = [1+1i, -1+1i, -1-1i, 1-1i] / sqrt(2);
M_qpsk = length(qpskConstellation);
k_qpsk = log2(M_qpsk);
qpskBits = [0 0; 0 1; 1 1; 1 0]; % Gray mapping

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
amplitudeTolerance = 0.25; %0,98, 0.8, 0.5, 0.25
phaseTolerance_deg = 3;  %30, 20, 5, 3

%% Colori simboli QPSK
myColor = [0.00, 0.45, 0.74;
           0.85, 0.33, 0.10;
           0.93, 0.69, 0.13;
           0.39, 0.83, 0.07];

%% Simulazione pattern e scatter plot

N_pattern = 0; %500;
tx_pattern = qpskConstellation(randi(4,N_pattern,1));

rx_ref_pattern = zeros(N_pattern,1);
rx_int_pattern = zeros(N_pattern,1);

W_used = zeros(ulaSize,N_pattern);
colors = zeros(N_pattern, 3);

% Prepara figura
fh = figure(1);
fh.WindowState = 'maximized';
sgtitle('Directional Modulation con ULA e QPSK', ...
        'Interpreter','Latex','FontSize',20);

subplot(2,2,1); hold on; grid on; box on; axis([-90 +90 -30 0]);
xlabel('$\varphi$ [deg]','Interpreter','Latex');
ylabel('Amplitude [dB]','Interpreter','Latex');

subplot(2,2,3); hold on; grid on; box on; axis([-90 +90 -180 +180]);
xlabel('$\varphi$ [deg]','Interpreter','Latex');
ylabel('Phase [deg]','Interpreter','Latex');

subplot(2,2,2); hold on; grid on; box on; axis square;
xlabel('I','Interpreter','Latex'); ylabel('Q','Interpreter','Latex');

subplot(2,2,4); hold on; grid on; box on; axis square;
xlabel('I','Interpreter','Latex'); ylabel('Q','Interpreter','Latex');

                                                                               
%% Vector generation
parfor n = 1:N_pattern
    [~,idx] = min(abs(tx_pattern(n) - qpskConstellation));
    colors(n,:) = myColor(idx,:);
    
    [w, ok] = FindSteeringVector(ulaSize, lambda, antennaLocation, ...
        theta_ref, phi_ref, azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg, B);
    
    fprintf("genero %d vettori \n", n);

    if ~ok
        warning('No solution found for the %d simbol \n', n);
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

subplot(2,2,4); hold on;
for n = 1:N_pattern
    plot(real(rx_int_pattern(n)), imag(rx_int_pattern(n)), 'o', ...
         'MarkerSize',6, 'MarkerFaceColor', colors(n,:), 'MarkerEdgeColor','k');
end

%% PLOT ampiezza e fase array
phi = linspace(-pi/2, pi/2, 1001);
theta = 0;

for n = 1:N_pattern
    [~, idxSym] = min(abs(tx_pattern(n) - qpskConstellation));
    thisColor   = myColor(idxSym,:);
    w = W_used(:,n);

    arrayResponse = zeros(1,length(phi));
    for idx = 1:length(phi)
        arrayResponse(idx) = (w'*exp(1i*myWaveVector(lambda,theta,phi(idx))' ...
            *antennaLocation).') / ulaSize;
    end

    % Ampiezza db
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
%   'FontSize',12, 'Callback', @(src,event) salvaVettori(W_used,ulaSize, B));


   % caso teorico
            % 3 antenne
            %w = [1 + 0i; 0.47618255771067436 - 0.8793464457948984i; -0.5465003434642405 - 0.8374588793448112i];  
            % 5 antenne 
            %w_theoretical = [1.0000 + 0.0000i; 0.4762 - 0.8793i; -0.5465 - 0.8375i; -0.9967 + 0.0818i; -0.4027 + 0.9153i];
            % 7 antenne
            %w_theoretical = [1.0000 + 0.0000i; 0.4762 - 0.8793i; -0.5465 - 0.8375i; -0.9967 + 0.0818i; -0.4027 + 0.9153i; 0.6132 + 0.7900i; 0.9866 - 0.1630i];
            

%% Simulazione BER
fileName = sprintf('steeringVector%d_%d-%d_%d_%d_%d.mat', ulaSize, grad, azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg, B);
if isfile(fileName)
    load(fileName, 'allVectors');
    nVectors = numel(allVectors);
    fprintf("Sono stati trovati %d set di vettori salvati per ULA=%d.\n", nVectors, ulaSize);
    %% Parametri BER
    EbN0_dB    = 0:1:9; %9     % intervallo Eb/N0
    nMax_error = 1000;      % numero massimo di errori prima di fermarsi
    
    % Risultati
    BER_dest_sim = zeros(size(EbN0_dB));
    BER_int_sim  = zeros(size(EbN0_dB));
    BER_int2_sim = zeros(size(EbN0_dB));

    %% Loop Eb/N0
    for k = 1:length(EbN0_dB)

        EbN0_lin = 10^(EbN0_dB(k)/10)

        N_symbols      = 0;
        N_errors_dest  = 0;
        N_errors_int   = 0;
        N_errors_int2  = 0;

        while N_errors_dest < nMax_error
            % Seleziona un vettore casuale dal file
            idxVec   = randi(nVectors);         
            W_test   = allVectors{idxVec};
            col      = randi(size(W_test,2));   
            w        = W_test(:,col);

            % steering vectors
            a_ref  = exp(1i*myWaveVector(lambda,theta_ref,phi_ref)'*antennaLocation).';
            a_int  = exp(1i*myWaveVector(lambda,theta_int,phi_int)'*antennaLocation).';
            a_int2 = exp(1i*myWaveVector(lambda,theta_int2,phi_int2)'*antennaLocation).';

            % random QPSK symbol
            symIdx = randi(4,1,1); 
            txSym  = qpskConstellation(symIdx);
            txBits = qpskBits(symIdx,:);


            gain_dest = abs((w'*a_ref)/ulaSize);
            gain_int  = abs((w'*a_int)/ulaSize);
            gain_int2 = abs((w'*a_int2)/ulaSize);
             
            rx_dest = gain_dest  * txSym;  
            rx_int  = gain_int   * txSym;
            rx_int2 = gain_int2  * txSym;
            

            % noise generation
            sigma2 = 1 / (EbN0_lin * k_qpsk); 
            
            noise_dest = sqrt(sigma2 / 2) * (randn + 1i*randn);
            noise_int  = sqrt(sigma2 / 2) * (randn + 1i*randn);
            noise_int2 = sqrt(sigma2 / 2) * (randn + 1i*randn);
            
            
            rx_dest_noise   = rx_dest + noise_dest;
            rx_int_noise    = rx_int  + noise_int;
            rx_int2_noise   = rx_int2 + noise_int2;
           
            % ML
            [~, idx] = min(abs(rx_dest_noise  - qpskConstellation));
            rxBits_dest = qpskBits(idx,:);

            [~, idx] = min(abs(rx_int_noise - qpskConstellation));
            rxBits_int = qpskBits(idx,:);

            [~, idx] = min(abs(rx_int2_noise  - qpskConstellation));
            rxBits_int2 = qpskBits(idx,:);

            % update counters variables
            N_symbols      = N_symbols + 1;
            N_errors_dest  = N_errors_dest + sum(txBits ~= rxBits_dest);
            N_errors_int   = N_errors_int  + sum(txBits ~= rxBits_int);
            N_errors_int2  = N_errors_int2 + sum(txBits ~= rxBits_int2);
        end

        % BER generation
        BER_dest_sim(k) = N_errors_dest / (k_qpsk*N_symbols);
        BER_int_sim(k)  = N_errors_int  / (k_qpsk*N_symbols);
        BER_int2_sim(k) = N_errors_int2 / (k_qpsk*N_symbols);

        fprintf("Eb/N0 = %d dB \t| Symbols: %d \t| BER_dest = %.4e\n", EbN0_dB(k), N_symbols, BER_dest_sim(k));
    end

    %% BER theory QPSK
    EbN0_lin    = 10.^(EbN0_dB/10);
    BER_theory  = qfunc(sqrt(2*EbN0_lin));
   
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
    %title(sprintf('BER QPSK, ULA\\_SIZE=%d', ulaSize), 'Interpreter', 'latex', 'FontSize',20);
    legend({'Theory','Intended receiver','Eavesdropper 1','Eavesdropper 2'},'Location','southwest', 'Interpreter','latex', 'FontSize', 16);
end




end


   

%% === Funzioni di supporto ===
function k = myWaveVector(lambda, theta, phi)
    k = 2*pi/lambda * [cos(theta)*cos(phi); cos(theta)*sin(phi); sin(theta)];
end
function [steeringVector, solutionIsFound] = FindSteeringVector(ulaSize, lambda, antennaLocation, ...
    theta_ref, phi_ref, azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg, B)

    maxIter = 1e8; 
    phi = linspace(-pi/2,pi/2,201); 
    theta = 0; 
    solutionIsFound = false;

    numLevels = 2^B;
    for n = 1:maxIter
        levels = (0:numLevels-1)'*(2*pi/numLevels);
        idx = randi(numLevels, ulaSize,1);
        phaseVector = levels(idx); % + deg2rad(2)*randn(ulaSize,1);
        steeringVector = exp(1i*phaseVector);
        
        %Theory 
        %phaseVector=2*pi*rand(ulaSize,1);
        %steeringVector=exp(1i*phaseVector); % random steering vector
       

        arrayResponse = zeros(1,length(phi));
        for index = 1:length(phi)
            arrayResponse(index) = (steeringVector' * ...
                exp(1i*myWaveVector(lambda,theta,phi(index))'*antennaLocation).') / ulaSize;
        end
        
        sel = (phi >= phi_ref-deg2rad(azimuthWidth_deg)/2) & ...
              (phi <= phi_ref+deg2rad(azimuthWidth_deg)/2);
        if ~any(sel), continue; end
        selectedArrayResponse = arrayResponse(sel);
        
        maxAmp_dB   = pow2db(abs(max(selectedArrayResponse)));
        minAmp_dB   = pow2db(abs(min(selectedArrayResponse)));
        angles_deg  = rad2deg(angle(selectedArrayResponse));
        
        if maxAmp_dB <= 0 && minAmp_dB >= -amplitudeTolerance && ...
           max(angles_deg) <= phaseTolerance_deg/2 && ...
           min(angles_deg) >= -phaseTolerance_deg/2
            solutionIsFound = true; 
            return;
        end
    end
end

function salvaVettori(W_used, ulaSize, B)
    fileName = sprintf('steeringVector%d_Quantization%d.mat', ulaSize, B);

    % Carica vettori esistenti se esistono
    if isfile(fileName)
        S = load(fileName, 'allVectors');
        if isfield(S,'allVectors')
            allVectors = S.allVectors;
            if ~iscell(allVectors)
                allVectors = {allVectors}; % forzato a cell array
            end
        else
            allVectors = {};
        end
    else
        allVectors = {};
    end

    % Salva ogni colonna di W_used come cella separata
    newVectors = num2cell(W_used, 1); % crea cell array riga
    newVectors = newVectors(:); % forza colonna

    % Forza allVectors ad essere cell array colonna
    if iscell(allVectors) && size(allVectors,2) > 1
        allVectors = allVectors(:);
    end

    % Concatenazione verticale
    allVectors = [allVectors; newVectors];

    save(fileName, 'allVectors', '-v7.3');
    fprintf('Vettori salvati in %s (totale set: %d)\n', fileName, numel(allVectors));
end
