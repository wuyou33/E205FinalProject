thetar = 1; % radians
thetardot = 0;
zeta = 1;
omega = 10;
Jhat = 1;
J = 2;
mup = 10;
mud = 10;

t = 0:0.01:10000;
[tout, ~, yout] = sim('satelliteProject', t);
goal = tout./tout.*thetar;

x1 = yout(:,5); 
x2 = yout(:,6); 
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
legend('Theta','Theta_m','Theta_r');
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


