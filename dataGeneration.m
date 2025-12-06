function test_ula_snr_v2_infinite
%% ===================== PARAMETRI =====================
close all; clc;

grad = -20;
theta_ref = deg2rad(0);
phi_ref   = deg2rad(grad);
azimuthWidth_deg   = 1;
amplitudeTolerance = 0.8;
phaseTolerance_deg = 20;
B = 2;
ulaSize = 7;

fileName = sprintf('steeringVector%d_%d-%d_%d_%d_%d.mat', ulaSize, grad,azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg, B);
N_pattern = 500;     % vettori generati per ciclo
nMax = 100000;       % limite massimo
pauseBetweenCycles = 1; % secondi di pausa tra cicli
fprintf("=== Generatore continuo ULA ===\nFile: %s | Max vettori: %d\n\n", fileName, nMax);

%% ===================== LOOP =====================
while true
   % Carica file esistente e NORMALIZZA allVectors come cell array colonna
   if isfile(fileName)
       S = load(fileName, 'allVectors');
       if isfield(S, 'allVectors')
           av = S.allVectors;
           % Se è numeric: split per colonne
           if isnumeric(av)
               allVectors = num2cell(av, 1);
               allVectors = allVectors(:); % colonna
           elseif iscell(av)
               % Garantisco che ogni cella contenga un singolo vettore colonna:
               tmp = {};
               for ii = 1:numel(av)
                   item = av{ii};
                   if isempty(item)
                       continue;
                   end
                   % se la cella contiene una matrice con più colonne,
                   % spezzala in tante celle (una per colonna)
                   if isnumeric(item) && size(item,2) > 1
                       cols = num2cell(item, 1);
                       cols = cols(:);
                       tmp = [tmp; cols]; %#ok<AGROW>
                   else
                       % altrimenti prendi l'elemento così com'è
                       tmp = [tmp; {item}]; %#ok<AGROW>
                   end
               end
               allVectors = tmp;
           else
               % Formato non previsto: inizia da zero
               allVectors = {};
           end
       else
           allVectors = {};
       end
   else
       allVectors = {};
   end
   % Assicurati che sia colonna
   if iscell(allVectors)
       allVectors = allVectors(:);
   else
       allVectors = {};
   end
   nCurrent = numel(allVectors);
   fprintf("🔹 Vettori nel file: %d\n", nCurrent);
   % Stop se raggiunto limite
   if nCurrent >= nMax
       fprintf("Raggiunto limite massimo (%d). Arresto.\n", nMax);
       break;
   end
   % Quanti ne possiamo ancora generare
   N_remaining = min(N_pattern, nMax - nCurrent);
   fprintf("▶ Generazione batch da %d vettori (al più)...\n", N_remaining);
   % Genera batch (parfor interno)
   W_used = generaVettoriSingoloCiclo(ulaSize, N_remaining, theta_ref, phi_ref, azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg, B);
   % Filtra solo quelli validi (non tutti i worker troveranno soluzione)
   validMask = any(W_used ~= 0, 1);
   W_valid = W_used(:, validMask);
   nValid = size(W_valid, 2);
   % Salva solo se ci sono vettori validi
   if nValid > 0
       newVectors = num2cell(W_valid, 1);
       newVectors = newVectors(:); % assicurati colonna
       % concatenazione verticale sicura
       allVectors = [allVectors; newVectors];
       save(fileName, 'allVectors', '-v7.3');
       fprintf("Salvati %d nuovi vettori validi (totale: %d)\n", nValid, numel(allVectors));
   else
       fprintf("Nessun vettore valido trovato in questo ciclo.\n");
   end
   pause(pauseBetweenCycles);
end
end
function W_used = generaVettoriSingoloCiclo(ulaSize, N_pattern, theta_ref, phi_ref, azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg, B)
f0 = 1575.42e6;
c  = 299792458;
lambda = c/f0;
spacing = lambda/2;
antennaLocation = zeros(3,ulaSize);
antennaLocation(2,:) = (0:ulaSize-1) * spacing;
W_used = zeros(ulaSize,N_pattern);
fprintf("  Calcolo %d vettori...\n", N_pattern);
tic;
parfor n = 1:N_pattern
   %fprintf("%d\n", n); 
   [w, ok] = FindSteeringVector(ulaSize, lambda, antennaLocation, ...
       theta_ref, phi_ref, azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg, B);
   if ok
       W_used(:,n) = w;
   end
end
fprintf("  Ciclo completato in %.1f s.\n", toc);
end


%% Funzioni di supporto
function k = myWaveVector(lambda, theta, phi)
   k = 2*pi/lambda * [cos(theta)*cos(phi); cos(theta)*sin(phi); sin(theta)];
end
function [steeringVector, solutionIsFound] = FindSteeringVector(ulaSize, lambda, antennaLocation, ...
   theta_ref, phi_ref, azimuthWidth_deg, amplitudeTolerance, phaseTolerance_deg, B)
   maxIter = 1e7;
   phi = linspace(-pi/2,pi/2,201);
   theta = 0;
   solutionIsFound = false;
   for n = 1:maxIter
       numLevels = 2^B;
       levels = (0:numLevels-1)'*(2*pi/numLevels);
       idx = randi(numLevels, ulaSize,1);
       phaseVector = levels(idx); %+ deg2rad(2)*randn(ulaSize,1);
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
       maxAmp_dB   = pow2db(abs(max(selectedArrayResponse)));
       minAmp_dB   = pow2db(abs(min(selectedArrayResponse)));
       angles_deg  = rad2deg(angle(selectedArrayResponse));
       if maxAmp_dB <= 0 && minAmp_dB >= -amplitudeTolerance && ...
          max(angles_deg) <= phaseTolerance_deg/2 && ...
          min(angles_deg) >= -phaseTolerance_deg/2
           steeringVector = steeringVector; %#ok<NASGU>
           solutionIsFound = true;
           steeringVector = steeringVector(:);
           return;
       end
   end
end
