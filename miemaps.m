function [EToN] = miemaps(Np, K)

% matrice con per ogni elemento il numero dei nodi 1 2 e 3
EToN = reshape(1:Np*K, Np, K)';

return;