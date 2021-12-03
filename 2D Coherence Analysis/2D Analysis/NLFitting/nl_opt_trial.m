close all; clear all; clc
T = linspace(80,2000,200);
P = multi_exp_with_oscillations_phi([80 300 1000],[1/50 1/100 1/200],T);
A = rand(2000,4)-0.5;
Y = zeros(length(T),size(A,1));
for k=1:size(A,1)
    Y(:,k) = P*A(k,:)'/range(P*A(k,:)') + rand(length(T),1)/10;
end
tic
[alpha, c, wresid, wresid_norm, y_est, Regression] =  ...
    varpro_multiway_2(Y,ones(length(T),size(A,1)),[50 200 1200]',4,@(x) multi_exp_phi(x,T),[0 200 800]',[250 900 2000]');
plot(Y); hold all; plot(y_est);
toc
figure;
 Z = Y-y_est;
 imagesc(abs(fft(Z,2^12,1)))