close;
Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 1000;             % Length of signal
t = (0:L-1)*T; 

A = zeros(768,1024);
B = zeros(768,1024);
C = zeros(768,1024);
D = zeros(1,1024);
for t = 1: 1024
    A(:,t)= (1+sin(pi*(t-256)^2/32768))/2;
    B(:,t)= (1+sin(pi*(t-512)^2/32768))/2;
    C(:,t)= (1+sin(pi*(t-768)^2/32768))/2;
    D(1,t)=sin(t^2);
end 
%{
figure
imshow(A);
figure
imshow(B);
figure
imshow(C);
figure
imshow((A+B+C)/3);
%}
Y = fft(A(1,:)+B(1,:)+C(1,:));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

