% clear
global phi A b N

% Define grid
h=0.1;
x=0:h:1;
y=0:h:1;

% Determine matrix dimensions
m=length(x);
N=m-2;
timer=tic;

% Boundary conditions
phi=ones(m);
phi(m,:)=2;

% Construct matrix A and vector b
A=-4*eye(N^2);
b=zeros(N^2,1);

for i=1:N
  for j=1:N
    row=(i-1)*N+j;

    add(row,i-1,j);
    add(row,i+1,j);
    add(row,i,j-1);
    add(row,i,j+1);
  end
end

% Solve for interior points
% X=inv(A)*b;
X=A\b;

% Combine boundary conditions and interior points
phi(2:m-1,2:m-1)=reshape(X,N,N);

% Plot result
mesh(x,y,phi)
loggedN(end+1)=N;
elapsed(end+1)=toc(timer);