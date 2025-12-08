alpha = 1/5; 
w0 = 2;

wr = w0/sqrt(1 - alpha^2);

e = exp(1);

phi = @(t) (e.^(-alpha*w0*t)).*sin(wr*t);
% t = 0:.01:50; 
p = 2500*pi/5103;
hold on
figure(1)
y = [];
for t = 0:.01:50
    if t <= 2*p 
        y = [y ,  0.1 * phi(t)]; 
    elseif t <= 4*p
        y = [y, 2 * phi(t)];
    else
        y = [y, phi(t)];
    end
end
t = 0:.01:50;
plot(t,y)
hold off

phi_tild = fft(y);

% G = @(w,w0)
