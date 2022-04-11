clear diff ddiff tdiff;

for i=1:5
    disp(eigenvalues(i))
end

for i=1:4
    diff(i)=eigenvalues(i)-eigenvalues(i+1)
end

for i=1:3
    ddiff(i)=diff(i)-diff(i+1)
end

for i=1:2
    tdiff(i)=ddiff(i)-ddiff(i+1)
end