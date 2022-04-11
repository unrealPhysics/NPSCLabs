clear;

i=1;

for e=-1:0.001:0
    A=shootingSchrodinger(e);
    E(i)=e;
    value(i)=A(1);
    i=i+1;
end

plot(E,value,'bo')

k=1;

for j=1:(i-2)
    if(value(j)*value(j+1)<0)
        eigenBounds=[E(j) E(j+1)];
        fullBounds(k,1)=eigenBounds(1);
        fullBounds(k,2)=eigenBounds(2);
        k=k+1;
    end
end