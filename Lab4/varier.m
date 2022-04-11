clear;

i=1;

for e=0:0.05:4
    A=shooting(e);
    E(i)=e;
    value(i)=A(1);
    i=i+1;
end

plot(E,value,'bo')