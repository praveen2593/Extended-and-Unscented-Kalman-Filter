% %(t,x,u,c) = needleModel(1,2,3,1/3);
% p1 = 1
% p2 = 1
% p3 = 1
% h1 = 2
% h2 = 2
% h3 = 2
% b1 = 3
% b2 = 3
% b3 = 4

u = [0.005;0;0;0]
c = 1/3

% Make up a trajectory save all points along the trajectory and
% measurements into a vector for use later.  This is because it is a test
% and we dont have an actual system to do this for us.
trajLength = 200;
X_IC = [2;3;4;0.8;0.6;0;0;0.8;-0.6];
trajectory = zeros(9,trajLength);
measurement = zeros(6,trajLength);
dT = 1;
IC = X_IC; % for ode looping
for i = 1:trajLength
    U = u;
    [~,y] = ode45(@(t,x) needleModel(t,x,U,c),[0;dT],IC);
    [P,Hs,B,Jx,Ju] = unpack(y(end,:));
    Xf = [P;Hs;B];
    
    trajectory(:,i) = Xf;
    measurement(1:3,i) = P + randn(3,1)*0.0005;
    measurement(4:6,i) = B + randn(3,1)*0.0001;
    IC = Xf;
end

% OK.  Now we have our test measurements...  So we need to iterate through
% all of the measuremnts using the EKF to estimate the state along the
% way...
filterOutput = zeros(9,trajLength);

X_l = X_IC + randn(9,1)*.1;
P_l = 0.01*eye(9);

for i= 1:trajLength
    U = u;
    Ym =  measurement(:,i);
    [X_o, P_o] = needleModelekf(Ym, U, X_l, P_l, dT);
    
    filterOutput(:,i) = X_o;
    
    X_l = X_o;
    P_l = P_o;
end

%%
figure(1)
subplot(2,1,1)
plot(trajectory(1,:),trajectory(2,:),'k',...
    measurement(1,:),measurement(2,:),'r.',...
    filterOutput(1,:),filterOutput(2,:),'b.');
title('XY trajectory');
subplot(2,1,2)
plot(dT*(0:(trajLength-1)),sqrt(sum((trajectory(1:3,:)-filterOutput(1:3,:)).^2)),'b.',...
    dT*(0:(trajLength-1)),sqrt(sum((trajectory(1:3,:)-measurement(1:3,:)).^2)),'r.')
title('XY Filtered Trajectory Error')
