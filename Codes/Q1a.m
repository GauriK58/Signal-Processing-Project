clear; close all; clc;

N = 200;             
Tsamp = 1;               
dT = Tsamp/10;        
Kmax = 4;                

n = 0:N-1;

sig1 = 0.9*sin(0.05*pi*n) + 0.4*cos(2*pi*n);      
sig2 = exp(-0.01*n) .* sin(0.1*pi*n);              
sig3 = 0.5*sin(0.12*pi*n)+0.3*sin(0.18*pi*n)+0.2*cos(0.06*pi*n);                       

k_n = randi([-Kmax Kmax], [1 N]);  
t_jit = n*Tsamp + k_n*dT;  

%jittered samples
jit1 = zeros(1, N);
jit2 = zeros(1, N);
jit3 = zeros(1, N);

for i = 1:N
    jit1(i) = sum( sig1 .* sinc( (t_jit(i)/Tsamp) - n ) );
    jit2(i) = sum( sig2 .* sinc( (t_jit(i)/Tsamp) - n ) );
    jit3(i) = sum( sig3 .* sinc( (t_jit(i)/Tsamp) - n ) );
end

%matrix
A = zeros(N, N);
for i = 1:N
    for m = 1:N
        A(i,m) = sinc( (t_jit(i)/Tsamp) - (m-1) );
    end
end

lambda = 0.001;

% A'Ax = A'y
% so x = (A'A)^-1 A'y
% x = (A'A + lambdaI)^-1 A'y cuz A'A might be singular

x_rec1 = (A'*A + lambda*eye(N)) \ (A'*jit1.');
x_rec2 = (A'*A + lambda*eye(N)) \ (A'*jit2.');
x_rec3 = (A'*A + lambda*eye(N)) \ (A'*jit3.');

%% MSE
MSE1 = mean((x_rec1.' - sig1).^2);
MSE2 = mean((x_rec2.' - sig2).^2);
MSE3 = mean((x_rec3.' - sig3).^2);

disp('MSE(signal 1) = '), disp(MSE1);
disp('MSE(signal 2) = '), disp(MSE2);
disp('MSE(signal 3) = '), disp(MSE3);

%% PLOTS
figure;

subplot(3,2,1)
plot(n, sig1); hold on;
plot(n, jit1, 'x');
title('Signal 1: True vs Jittered Samples');
legend('True','Jittered');

subplot(3,2,2)
stem(n, x_rec1);
title('Signal 1: Reconstructed');

subplot(3,2,3)
plot(n, sig2); hold on;
plot(n, jit2, 'x');
title('Signal 2: True vs Jittered Samples');
legend('True','Jittered');

subplot(3,2,4)
stem(n, x_rec2);
title('Signal 2: Reconstructed');

subplot(3,2,5)
plot(n, sig3); hold on;
plot(n, jit3, 'x');
title('Signal 3: True vs Jittered Samples');
legend('True','Jittered');

subplot(3,2,6)
stem(n, x_rec3);
title('Signal 3: Reconstructed');
