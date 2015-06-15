y = 3*ones(3001,1) + 2*(randn(3001,1)); % signal amplitude is 3 units and noise amplitude is 2 unit
y_p(1) = 0;
k = 0.5.*(2*rand(403,1)-1);
psi(:,1) = exp(-(((-200-1:200+1)-30).^2)/(2*1000))*(1/sqrt(2*pi*1000));
plot((-200-1:200+1)*0.1,psi(:,1))
%%
for t = 1:3000
    t
    for ep = 1:2
        V = -(480*(y(t) - y_p(t))).*k;
        psi(2:403-1,t+1) = psi(2:403-1,t) + ((1i*0.01)/(2*0.1*0.1*0.1)).*(psi(3:403,t) - 2*psi(2:403-1,t) + psi(1:403-2,t)) - 1i.*0.01.*V(2:403-1).*psi(2:403-1,t);
        psi(:,t+1) =  psi(:,t+1)./sum(abs(psi(:,t+1)));
        p(:,t+1) = abs(psi(:,t+1)).^2/(sum(abs(psi(:,t+1)).^2));
        y_p(t+1) = (-200-1:200+1)*0.1'*p(:,t+1);
        
        k(:) = k(:) + (0.01*0.08*(y(t) - y_p(t))).*p(:,t+1);
        
    end
    plot((-200-1:200+1)*0.1,smooth(smooth(abs(psi(:,t+1)))))
    pause(0.001)
end
figure(2)
plot([1:3001],y,'-b',[1:3001],y_p,'-r',[1:3001],3*ones(3001,1),'-g')
figure(4)
plot((-200-1:200+1)*0.1,p(:,3001))
