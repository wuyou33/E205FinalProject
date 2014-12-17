clear 

zeta = 1;
omega = 10;
Jhat = 1;
J = 2;
mup = 100;
mud = 100;
Amplitude = 1;
Frequency = 1;
ditherAmp = 0.00;
ditherFreq = 100;
Period = 5;
inputChoice = 2;

initialAngle = 0;
initialVelocity = 0;

t = 0:0.01:100;
[tout, ~, yout] = sim('satelliteProject', t);
% goal = tout./tout.*thetar;

x1 = yout(:,5); 
x2 = yout(:,6); 
x3 = yout(:,9); 
x4 = yout(:,10); 
x3dot = yout(:,11); 
x4dot = yout(:,12); 
x2dot = yout(:,13); 
e = yout(:,14); 
edot = yout(:,15); 

u = yout(:,7); 
input = yout(:,8); 

figure(2)
clf
plot(tout,input);

T = 0.01; %s
I = T*sum(x1.^2+x2.^2+(u.^2)/100);

disp('I: ')
disp(I);

figure(1)
clf
subplot(3,1,1);
plot(tout,yout(:,1)); % Theta
hold all
plot(tout,yout(:,4)); % Theta_m
% plot(tout,goal,'--'); % Theta_r
% ylim([0.04 0.07]);
xlabel('Time (s)');
ylabel('Theta (rad)');
legend('Theta','Theta_m');
title('Theta over time');
subplot(3,1,2);
plot(tout,yout(:,2)); % P
hold all
% plot(tout,tout./tout*J*omega^2); % P_hat
xlabel('Time (s)');
ylabel('P');
title('P over time');
subplot(3,1,3); % D
plot(tout,yout(:,3));
hold all
% plot(tout,tout./tout*J*2*zeta*omega); % D_hat
xlabel('Time (s)');
ylabel('D');
title('D over time');
% subplot(4,1,4); % u
% plot(tout,u);
% xlabel('Time (s)');
% ylabel('u(t)');
% title('u over time');

a1 = omega^2;
ap = 1/(J*mup);
ad = 1/(J*mud);
V = a1*x1.^2+x2.^2+ap*x3.^2+ad*x4.^2;

figure(4)
clf
hold all
plot(tout, V);
Vdot = -4*zeta*omega*x2.^2;
plot(tout, Vdot);
% Vdot2 = 2*a1*x1.*x2+2*x2.*x2dot+2*ap*x3.*x3dot+ad*x4.*x4dot;
% plot(tout, Vdot2)
% x1zero = (a1-omega^2)*2.*x1.*x2;
% x2nonzero = -4*zeta*omega*x2.^2;
% x3zero = (x2.*e./J + ap*x3dot).*2.*x3;
% x4zero = x2.*edot./J + ad*x4dot.*2.*x4;
% plot(tout, x1zero+x2nonzero+x3zero+x4zero);

% plot(tout, a1*x1.^2);
% plot(tout, x2.^2);
% plot(tout, ap*x3.^2);
% plot(tout, ad*x4.^2);
legend('v','vdot');
% ylim([-5 5]);
%, 'vdot2','vdot3','1','2','3','4'

% % x1zero = (a1-omega^2)*2.*x1.*x2;
% % x2nonzero = -4*zeta*omega*x2.^2;
% % x3zero = (x2.*e./J + ap*x3dot).*2.*x3;
% % x4zero = x2.*edot./J + ad*x4dot.*2.*x4;
% figure(5)
% clf
% plot(tout, x1zero)
% hold all
% plot(tout, x3zero)
% plot(tout, x4zero)
% legend('x1','x3', 'x4');
% 
% 
% x2dotBright = e.*x3/J+edot.*x4/J-2*zeta*omega.*x2-omega^2.*x1;
% figure(6)
% clf
% plot(tout, 4*x2dot);
% hold all
% plot(tout,x2dotBright);
% plot(tout, x2);
% plot(tout, x1);
% ylim([-20 20]);
% legend('x2dot','x2dotBright','x2','x1');
% 
% % figure(7)
% % clf
% % plot(tout, e)
% % hold all
% % plot(tout, edot)
% % legend('e','edot');