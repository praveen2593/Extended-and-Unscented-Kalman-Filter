function [X_o, P_o] = needleModelekf(Ym, U, X_l, P_l, dT)
c = 1/3; % curvature coefficient
Q = 0.00001*eye(9); %Process noise covariance
R = diag([.0005 .0005 .0005 .0001 .0001 .0001]); %Measurement noise covariance
H_o = [eye(3) zeros(3,6); zeros(3,6) eye(3)];


% for all of our measurments
IC = [X_l;reshape(eye(9),81,1);zeros(36,1)];
[~,y] = ode45(@(t,x) needleModel(t,x,U,c),[0;dT],IC,odeset('AbsTol',1e-9));
[P,Hs,B,Jx,~] = unpack(y(end,:));
Xf = [P;Hs;B];  
I = eye(9);

Ye = H_o*Xf; %X is pulled from "repack.m"

P_x = Jx*P_l*Jx'+Q;  %P_x = df_dx*P*df_dx'+Q %State covariance
K = P_x*H_o'*(H_o*P_x*H_o'+R)^-1; %Kalman Gain. Based on H_o and Z, should H_o here become Z?
P_o = (I - K*H_o)*P_x; %Output covariance. Should H_o become Z?
%g(X_next) = X %Measurement function?
X_o = Xf+K*(Ym - Ye); %Output state. Where Z is current measurement, 'g' is the measurement function. Parenthases should be measured value minus predicted value. Initial equation had X_o = Xf+K*(Z - Ym) but Z & Ym are equal hence noise had to be added to measurement

end


% save X_e as we go

%end for

% now plot actuals and expected to see how we did...

