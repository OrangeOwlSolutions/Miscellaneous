clear all
close all
clc

load Source_Matrix.txt
load U.txt
load V.txt
load S.txt

Source_Matrix = reshape(Source_Matrix, 8, 8);
U             = reshape(U, 8, 8);
V             = reshape(V, 8, 8);

[UU, SS, VV]  = svd(Source_Matrix);

100 * sqrt(sum(abs(diag(SS) - S).^2) / sum(abs(diag(SS)).^2))

100 * sqrt(sum(sum((abs(Source_Matrix - UU * SS * VV.').^2))) / sum(sum((abs(Source_Matrix).^2))))

100 * sqrt(sum(sum((abs(Source_Matrix - U  * diag(S)  * V).^2))) / sum(sum((abs(Source_Matrix).^2))))

