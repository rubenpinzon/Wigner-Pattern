function W = filterbins(W,start, finish)
%
%       Auxiliar function to extract bins in the interval [start, finish]
%
%Ruben Pinzon

len_w = length(W);
for n = 1 : len_w
   W(n).y  = W(n).y(:,start+1:finish); 
   W(n).T  = finish-start;
   W(n).binsUsed = [start+1, finish]; 
end