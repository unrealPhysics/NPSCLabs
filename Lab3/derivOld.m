function value=derivOld(x)
global g L

value(1)=x(2);
value(2)=-g/L * sin(x(1));