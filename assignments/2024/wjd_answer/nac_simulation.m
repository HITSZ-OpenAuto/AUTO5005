%% params
T=0.001;N=10000;
K1=eye(2,2);K2=eye(2,2);
t=0:T:T*(N-1);
m1=1;m2=1.5;
l1=0.2;l2=0.3;lc1=0.1;lc2=0.15;
J1=0.013;J2=0.045;
g=9.8;
theta=[m1*lc1*lc1+m2*(l1*l1+lc2*lc2)+J1+J2 m2*l1*lc2 m2*lc2*lc2+J2 m1*lc1+m2*l1 m2*lc2];

%% goals
I=ones(1,N);
% case 1
% qd=[0.3*I; 0.1*I];dqd=zeros(2,N);ddqd=zeros(2,N);
% case 2
qd=[(0.3+0.02*sin(t)); (0.1+0.01*cos(t))];dqd=[0.02*cos(t);-0.01*sin(t)];
ddqd=[-0.02*sin(t);-0.01*cos(t)];
q=zeros(2,N);dq=zeros(2,N);ddq=zeros(2,N);
%% calculation
for i=1:N-1
M=[theta(1)+2*theta(2)*cos(q(2,i)) theta(3)+theta(2)*cos(q(2,i));theta(3)+theta(2)*cos(q(2,i)) theta(3)];
C=[-theta(2)*sin(q(2,i))*dq(2,i) -theta(2)*sin(q(2,i))*(dq(2,i)+dq(1,i));...
    theta(2)*sin(q(2,i))*dq(1,i) 0];
ddq(:,i+1)=inv(M)*(-K1*(q(:,i)-qd(:,i))-C*dq(:,i)-K2*(dq(:,i)-dqd(:,i))+C*dqd(:,i)+M*ddqd(:,i));
% note that there's no g in the expression, because it's eliminated
q(:,i+1)=dq(:,i)*T+q(:,i);
dq(:,i+1)=ddq(:,i)*T+dq(:,i);
end
e=q-qd;de=dq-dqd;
figure;
subplot(1,2,1);plot(t,e(1,:),t,e(2,:));
legend('error of $q_1$','error of $q_2$','Interpreter','latex');
subplot(1,2,2);plot(t,de(1,:),t,de(2,:));
legend('error of $\dot{q}_1$','error of $\dot{q}_2$','Interpreter','latex');
% sgtitle('$q_d=[0.3,0.1],\dot{q}_d=[0,0],q(0)=\dot{q}(0)=[0,0]$','Interpreter','latex')
sgtitle('$q_d=[0.3+0.02\sin t,0.1+0.01\cos t],\dot{q}_d=[0.02 \cos t,-0.01\sin t],q(0)=\dot{q}(0)=[0,0]$','Interpreter','latex')