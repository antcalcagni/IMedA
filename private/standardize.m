function [XC] = standardize(X)

% IN :
% X Matrice da Centrare e Standardizzare (n individui x m variabili)
%
% OUT : 
% Matrice XC, dove XC'*XC=I


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,p]=size(X);                      %DIMENSIONI DELLA MATRICE n=Righe; p=Colonne        %
PS=eye(n)/sqrt(n-1);                %PESI PARI A 1/n-1                                  %
M=(diag(std(X)))^-1;                %METRICA PARI A 1/S                                 %
C=X-ones(n,p)*diag(mean(X));        %CENTRA                                             %
%XC=PS*C*M;                          %MATRICE CENTRATA E STANDARDIZZATA                  %
XC=C*M;                          %MATRICE CENTRATA E STANDARDIZZATA                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
