function add(row,x,y)
global phi A b N

if min(x,y)<1 | max(x,y)>N
   b(row)=b(row) - phi(y+1,x+1);
else
   A(row,(x-1)*N+y)=1;
end

