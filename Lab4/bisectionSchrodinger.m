function midpoint = bisectionSchrodinger(lower,upper)
format compact

for i=1:20
  midpoint=0.5*(lower+upper);
  disp([i midpoint])

  if shootingSchrodinger(midpoint)*shootingSchrodinger(upper)<0
     lower=midpoint;
  else
     upper=midpoint;
  end
end
