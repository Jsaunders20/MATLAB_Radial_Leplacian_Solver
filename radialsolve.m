% Consistency check
% A manufactured solution u = x^2

% boundary data
G = R.^2.*cos(THETA).^2; 
% L u = f
F = 2 + 0*R;                      

RHS = G.*(DBDY ~= 0) + F.*(DBDY == 0);

rhs = reshape(RHS,N,1);

u =A\rhs; %u=A^-1 * rhs

U = reshape(u,size(THETA));
U_exact = R.^2.*cos(THETA).^2;


% Display the surface plot in Cartesian coordinates
surf([X;X(1,:)], [Y;Y(1,:)], [U;U(1,:)]- [U_exact;U_exact(1,:)], 'edgecolor', 'none');
title('Polar Surface Plot');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');

axis equal
