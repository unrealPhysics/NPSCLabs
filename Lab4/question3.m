clear;

varierSchrodinger
for k=1:3
    eigenvalues(k)=bisectionSchrodinger(fullBounds(k,1),fullBounds(k,2))
end
% % question3plotter
% for k=1:7
%     fraction(k)=eigenvalues(k)/eigenvalues(k+1)
% end