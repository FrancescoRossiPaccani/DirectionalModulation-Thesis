% Conta il numero di vettori unici nel file allVectors

clear; clc;

fileName = 'steeringVector3_-20-1_8.000000e-01_20_6.mat';
tol = 0;  % tolleranza numerica

fprintf("=== Conteggio vettori unici in '%s' ===\n", fileName);

if ~isfile(fileName)
    error('File non trovato.');
end

S = load(fileName, 'allVectors');
if ~isfield(S, 'allVectors')
    error("Il file non contiene la variabile 'allVectors'.");
end

allVectors = S.allVectors;
if ~iscell(allVectors)
    error("'allVectors' non è una cell array.");
end

N = numel(allVectors);
fprintf("Totale vettori nel file: %d\n", N);

% --- Conversione in matrice ---
vecList = cellfun(@(v) v(:).', allVectors, 'UniformOutput', false);
V = cell2mat(vecList);

% --- Confronto ---
Vcat = [real(V), imag(V)];
[~, ia] = uniquetol(Vcat, tol, 'ByRows', true);

nUnique = numel(ia);

fprintf("Numero di vettori unici (entro tol = %.1e): %d\n", tol, nUnique);
