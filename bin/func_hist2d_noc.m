function [Matrix_Counts] = func_hist2d(DataX,DataY,BinsX,BinsY);
%Histogramme 2D
%Input:
% Les données, sous forme de couple DataX (abscisses) et DataY (ordonnées)
% Les bins, sous forme de vecteurs avec BinsX = vecteur contenant les edges d'abscisses
% et  BinsY = vecteur contenant les edges d'ordonnées
%Output:
% Matrix_Counts, matrice 2D de dimension BinsX * BinsY, qui compte le nombre d'éléments 
% (DataX,DataY) qui tombent dans chaque case (BinX,BinY).
%P. Bottelin, Oct 2012; modified from Doug Hull, Matlab

%Mise en vecteurs colonne
% DataX=DataX(:);
% DataY=DataY(:);
% BinsX=BinsX(:);
% BinsY=BinsY(:);

NearestX=interp1(BinsX,1:numel(BinsX),DataX,'nearest');
NearestY=interp1(BinsY,1:numel(BinsY),DataY,'nearest');


Matrix_Counts=(accumarray([NearestX NearestY],1,[length(BinsX) length(BinsY)]))';
end