clear all
%EXAM 2
%Question 3

Lx = 15;
x = linspace(0,Lx,30);
Nx = length(x);
dt = 0.1;
dx = x(2) - x(1);
% a = 1; %alpha
tn = 100;

% T_steady = x.^2 .* exp(-x);
S_x = -(x.^2 - 4.*x + 2) .* exp(-x);

%INITIAL AND BOUNDARY CONDITIONS

T = zeros(Nx, 1);
TT = zeros(Nx, 1); 

TT(1,:) = 0;
TT(Nx,:) = (Lx^2)*exp(-Lx);

%EXPLICIT EULER

for iterations = 1:tn
    for i = 2:Nx-1;
        T(i) = 0.5*(TT(i+1) + TT(i-1)) + dx^2*0.5*(S_x(i));
    end
    TT = T;
end

%CRANK-NICHOLSON
% 
% A_xU(i) = (T(i+1) - 2*T(i) + T(i-1))/dx^2;
% 
% for time = 0:dt:20
%     for i = 2:Nx-1;
%         TTT = TT(i) - ((1/2) * (A_xU(i+1) + A_xU(i)) * dt);
%     end
% end
%     

%EXACT SOLUTION

for i = 1:Nx;
    T_exact(i) = x(i).^2 * exp(-x(i));
end

figure(1);
plot(T_exact, 'LineWidth', 2)
hold on
plot(T, 'LineWidth', 2)
title('Nx = 10')
legend('Exact', 'Explicit Euler')
hold off
