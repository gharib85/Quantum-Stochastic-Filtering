clc
close all
clear all

del_x=0.1;
del_t=0.01;
t=0:del_t:10;
[r,epoch]=size(t);
N=400;
h=1;
K=rand(N,epoch);
v=zeros(N,epoch);
shi=rand(N,epoch);
shi(:,1)=shi(:,1)/norm(shi(:,1)); 
beta=0.08;
m=0.1;
gamma=2;
xi=480;
yhat=zeros(1,epoch);

% Signal
x=2;
y=awgn(x,6);

a=zeros(N,1);
b=zeros(N,1);
nu=y-yhat(1);
% Learning
for i=1:epoch-1
    
   for I=1:gamma
       i
       v(:,i)=-(xi*nu*K(:,i));
       for k=2:N-1
          a(k)=-del_t*v(k,i)*shi(k,i) + ((del_t/(2*m*del_x*del_x))*(shi(k+1,i)-2*shi(k,i)+shi(k-1,i)));
       end
       a(1)=-del_t*v(1,i)*shi(1,i) + ((del_t/(2*m*del_x*del_x))*(shi(2,i)-2*shi(1,i)));
       a(N)=-del_t*v(N,i)*shi(N,i) + ((del_t/(2*m*del_x*del_x))*(-2*shi(N,i)+shi(N-1,i)));
       for k=1:N
          b(k)=a(k)*1j;
       end
       for k=1:N
          shi(k,i+1)=shi(k,i)+b(k); 
       end
       shi(:,i+1)=(shi(:,i+1)/(norm(shi(:,i+1))));
       X=-199:1:200;
       s=0;
       for k=1:N
          s=s+(X(k)*del_x*((norm(shi(k,i+1)))^2)); 
       end
       yhat(i+1)=s;
       nu=y-yhat(i+1);
       for k=1:N
          K(k,i+1) = K(k,i)+beta*nu*((norm(shi(k,i+1)))^2); 
       end
   end
end


% Plotting training signals
figure(1)
plot(t,x,'b');
hold on;
plot(t,y,'r');
hold on;
plot(t,yhat,'g');

% Making pdf of training signal
pd=makedist('Normal');
y=-20:0.01:20;
pdf_y=pdf(pd,y);
x=-20:0.01:20;
pdf_x=pdf(pd,x);

%Plotting pdf of training signal
figure(2)
plot(y,pdf_y,'b');
hold on;
plot(x,pdf_x,'r');