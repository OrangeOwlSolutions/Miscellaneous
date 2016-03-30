clear all
close all

Nrows           = 10;
Ncols           = 8;
numMatrices     = 1024;

A           = load('inputMatrices.txt');
S           = full_svd(A, 10, 8, numMatrices);
Sapproach   = load('singularValues.txt');

RelativeErrorPercentage = 100 * (sum((S - Sapproach).^2) / sum(S.^2))
