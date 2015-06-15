%%
clc
clear all
close all


n_d = 200;
n_limit = 200; 
n_n = 2*n_limit+3;

y = 5*ones(n_d,1) + (0*randn(n_d,1));
y_p = zeros(n_d,1);
k = 0.5.*(2*rand(n_n,1)-1);
V = zeros(n_n,1);
psi = zeros(n_n,n_d+1) + 1i*zeros(n_n,n_d+1);
p = zeros(n_n,1) ;
del_x = 0.1;
del_t = 0.01;
x = (-n_limit-1:n_limit+1)*del_x;
tau = 480;
m = 0.1;
B = 0.07;
m_ep = 5;

%save('data1.mat')

%%


%load('data1.mat')
for s = 1:n_n
    psi(s,1) = exp(-((s-5)^2)/(2*200))*(1/sqrt(2*pi*200));
end
figure(1)
plot(x,psi(:,1))
grid on
title('initial packet')


%%


for t = 1:n_d
    fprintf('data %d/%d \n',t,n_d)
    y_p(t) = 0;
    for ep = 1:m_ep
        err = (y(t) - y_p(t));
        V = -(tau*err).*k;
    subplot(2,2,1)
    plot(x,k)
    title('k')
    subplot(2,2,2)
    plot(x,V)
    title('v')
        for l = 2:n_n-1
            psi(l,t+1) = psi(l,t) + ...
            ((1i*del_t)/(2*m*del_x*del_x)).*(psi(l+1,t) - 2*psi(l,t) + psi(l-1,t))...
            - 1i*del_t*V(l)*psi(l,t);
        end
        
        psi(:,t+1) =  psi(:,t+1)./sum(abs(psi(:,t+1)));
        
        psi(1,t+1) = psi(2,t+1);
        psi(n_n,t+1) = psi(n_n-1,t+1);


        p(:,t+1) = abs(psi(:,t+1)).^2;
        p(:,t+1) = smooth(p(:,t+1));
       subplot(2,2,3)
        plot(x,p(:,t+1))
        %psi(:,t+1) = p(:,1);
        p(:,t+1) = p(:,t+1)./(sum(p(:,t+1)) + 0.000001);
        %[mx(t),id(t)] = max(p(:,t+1));
        y_p(t) = 0;
                    
        for r = 1:n_n
           y_p(t) = y_p(t) + x(r)*p(r,t+1);
        
           k(r) = k(r) + del_t*B*err*p(r,t+1);
        end
    end
    
    %plot(x,p(:,t+1),'-')


    subplot(2,2,4)
    plot(y_p,'-')
    title(t)

    pause(0.001);

end

%%
title('p')
figure(30)
plot(y_p)

title('yp')
figure(400)
plot(y,'-b')
hold on
plot(y_p,'-r')
plot(5*ones(n_d,1),'-g')
hold off

title('y')

sum(p(:,n_d+1))