clear;

N = 100;

theta = linspace(0,2*pi,N+1);
deltr = 2*pi/N;
theta = theta(1:end-1);

RADIUS = 6;
innerradius = 2;

% r = innerradius:deltr:RADIUS;
r= sin(2.*theta):deltr:RADIUS;

%capital for gridded, lower for _____
[THETA,R] = ndgrid(theta,r);   

% anytime you have a componenet-wise/element wise calculation you need to
% have a . EX: .^
%surf(X,Y,X.^2+Y.^2);

Ntheta = length(theta); 
Nr = length(r);

i = 1:Ntheta;
j = 1:Nr;

[I,J] = ndgrid(i,j);

N = Ntheta*Nr; %total number of points on the gird

index = reshape(1:N,size(I));

A = sparse(N,N); %N*N matrix of 0's

% RJR interior vertices
for i = 2:Ntheta-1
    j = 2:Nr-1;
    center = index(i,j);
    out    = index(i,j+1);
    in     = index(i,j-1);
    left   = index(i-1,j);
    right  = index(i+1,j);

    rj = r(j); % Radial positions for the current set of points

    % Centered difference for polar Laplacian
    %
    % u_rr + u_r/r + u_tt/r^2 = (u(i,j+1)+u(i,j-1)-2*u(i,j))/dr^2
    %                         + (u(i,j+1)-u(i,j-1))/dr/r
    %                         + (u(i+1,j)+u(i-1,j)-2*u(i,j))/dt^2/r^2
    %                         = u(i,j+1)*(1/dr^2 + 1/dr/r)
    %                         + u(i,j-1)*(1/dr^2 - 1/dr/r)
    %                         + u(i+1,j)/dt^2/r^2
    %                         + u(i-1,j)/dt^2/r^2
    %                         + u(i,j)*(-2/dr^2-2/dt^2/r^2)

    A(sub2ind(size(A), center, center)) = (-2./deltr^2 - 2./(rj.^2 * deltr^2));
    A(sub2ind(size(A), center, out))    = 1./deltr^2 + 1./(2*deltr*rj);
    A(sub2ind(size(A), center, in))     = 1./deltr^2 - 1./(2*deltr*rj);
    A(sub2ind(size(A), center, left))   = 1./(rj.^2 * deltr^2);
    A(sub2ind(size(A), center, right))  = 1./(rj.^2 * deltr^2);

end

for j = 2:Nr-1
    i = 1;
    center = index(i,j);
    out    = index(i,j+1);
    in     = index(i,j-1);
    left   = index(Ntheta,j);
    right  = index(i+1,j); 
    
    rj = r(j); 

    A(center, center) = (-2./deltr^2 - 2./(rj.^2 * deltr^2));
    A(center, out)    = 1./deltr^2 + 1./(2*deltr*rj);
    A(center, in)     = 1./deltr^2 - 1./(2*deltr*rj);
    A(center, left)   = 1./(rj.^2 * deltr^2);
    A(center, right)  = 1./(rj.^2 * deltr^2);

end

for j = 2:Nr-1
    i = Ntheta;
    center = index(i,j);
    out    = index(i,j+1);
    in     = index(i,j-1);
    left   = index(i-1,j);
    right  = index(1,j); 
    
    rj = r(j); 

    A(center, center) = (-2./deltr^2 - 2./(rj.^2 * deltr^2));
    A(center, out)    = 1./deltr^2 + 1./(2*deltr*rj);
    A(center, in)     = 1./deltr^2 - 1./(2*deltr*rj);
    A(center, left)   = 1./(rj.^2 * deltr^2);
    A(center, right)  = 1./(rj.^2 * deltr^2);

end

B=A;
amp = 2;
period = 3;
DBDY = zeros(Ntheta, Nr);    % Initialize DBDY as a matrix of zeros (double)
tolerance = 5 * deltr;

condition = find(R <=  innerradius + abs(amp*sin(period.*THETA))); 
DBDY(condition) = 1;

condition = find(R == max(max(R)));
DBDY(condition) = 2 ;

[X, Y] = pol2cart(THETA, R);
surf(X,Y, DBDY);

title('Boundary Condition Application');
xlabel('X');
ylabel('Y');
zlabel('Boundary Applied');

% Set rows of A corresponding to a boundary element as 1 or 0
ind = find(DBDY);
for i = 1:length(ind)
    k = ind(i);
    A(k,:) = 0;
    A(k,k) = 1;
end









