%Note: Code behaves eratic when the value of numt, specified in the code is
%high. Make sure 'Clear all' command is executed before running the
%code.Clear all has not been hardcoded because it decreases the timing
%performance of the code. Hence, when necessary.

%% System Inputs initialization
v=3;
omega=[1;2;3];
%% Parameters Initialization
numt = 150; % time steps of the trajectory
tspan=[0 0.025]; % time period of integration
d=9; % Dimensionality of the state vector
muBar_state=zeros(9,1,numt); % Mean state Vector
real_obv=zeros(6,numt); % Measurements from the real system
state_trajec=zeros(9,numt); % State Measurements from the real system
observation=zeros(6,1,2*d+1,numt); % Observation vectors for each sigma point
sv_sp=zeros(9,1,2*d+1); % State vectors for each sigma point
sp=zeros(9,1,2*d+1); % Sigma Points
sp_tr=zeros(9,1,2*d+1); % Sigma Points after the first transformation
muBar_obv=zeros(6,1,numt); % Mean observation vector
InoCov=zeros(6,6); % Innovation Covariance Martrix
SO_Cov=zeros(9,6); % State-Observation Covariance matrix
mu_t=zeros(9,numt); % Mean State vector matrix
a= zeros(9,1,1,1); % Buffer variable for Mean of State vector Sigma Points
b= zeros(9,9);     % Buffer variable for Covariance of the transformed Sigma Points
q= zeros(6,1,1,1); % Buffer variable for Mean Observation
L=zeros(6,6,1,1);  % Buffer variable for Innovation covarianc
x1=zeros(9,6);     % Buffer variable for State-Observation Covariance
R=0.01*eye(9); % Noise covariance of process
Q=diag([.001,.001,.001,.001,.001,.001]); % Noise covariance of observation (measurement)
mu_treal=[1;3;4;0;1;0;0;1;0]; % Actual initial state vector

%%  trajectory Generator
for trajec=1:numt %%iteration for each sigma point
    if (trajec==1)
        state_trajec(1:9,trajec)=mu_treal(:,1); %%state vector at time =1
        
    else
        
        [jk,ele]=ode45(@(trajec,x)statevecInteg(trajec,x,v,omega),tspan,state_trajec(1:9,trajec-1));
        state_trajec(1:9,trajec)=ele(end,:)';
        
    end
    
    real_obv(:,trajec)=observe(state_trajec(1:3,trajec),state_trajec(7:9,trajec)) + Q*randn(6,1);
end
%% weights and Covariance Initialization
[wt_pt,wt_cov]=sigmaPntCovWt(d);
covariance(1:9,1:9)= eye(9);
mu_t(1:9,1) = state_trajec(1:9,1)+randn(9,1);
%% Filter
for t=2:numt %% iterations for time steps equal to 0.025 seconds
    %% Sigma Points for each state vector
    sp(1:9,1,1:2*d+1)=sigmapnt(mu_t(1:9,t-1),covariance(1:9,1:9),d);
    
    %% Vector rates integration
    
    for sgmp=1:2*d+1 %%iteration for each sigma point
        [xy,element]=ode45(@(t,x)statevecInteg(t,x,v,omega),tspan,sp(1:9,1,sgmp));
        sv_sp(1:9,1,sgmp)=element(end,:)';
    end
    
    %% Mean of State vector Sigma Points
    a = zeros(9,1,1,1);
    for msgmp=1:2*d+1
        a = a+wt_pt(msgmp)*sv_sp(1:9,1,msgmp);
    end
    muBar_state(:,1,t) = a;
    
    %% Covariance of the transformed Sigma Points
    b= zeros(9,9);
    for msgmcov=1:2*d+1
        covariance(1:9,1:9)= b+wt_cov(msgmcov)*(sv_sp(1:9,1,msgmcov)-muBar_state(1:9,1,t))*(sv_sp(1:9,1,msgmcov)-muBar_state(1:9,1,t))';
        b=covariance(1:9,1:9);
    end
    covariance(:,:)=covariance(:,:)+R;
    %% Sigma points around the transformed mean
    
    sp_tr(1:9,1,1:2*d+1)=sigmapnt(muBar_state(1:9,1,t),covariance(1:9,1:9),d);
    
    
    %% Observation matrix determination
    for sgmp_tr=1:2*d+1
        observation(:,1,sgmp_tr,t)=observe(sp_tr(1:3,1,sgmp_tr),sp_tr(7:9,1,sgmp_tr));
    end
    %% Mean Observation
    q=zeros(6,1,1,1);
    for msgm_obv=1:2*d+1
        q=q+wt_pt(msgm_obv)*observation(1:6,1,msgm_obv,t);
    end
    muBar_obv(:,1,t)=q;
    %% Innovation covariance
    L=zeros(6,6,1,1);
    for f=1:2*d+1
        L=L+ wt_cov(f)*(observation(:,1,f,t)-muBar_obv(:,1,t))*(observation(:,1,f,t)-muBar_obv(:,1,t))';
    end
    InoCov(:,:)=L;
    InoCov(:,:)=InoCov(:,:)+Q ;
    
    %% State-Observation Covariance
    x1=zeros(9,6);
    for g=1:2*d+1
        x1= x1+ wt_cov(g)*(sp_tr(:,1,g)-muBar_state(:,1,t))*(observation(:,1,g,t)-muBar_obv(:,1,t))';
    end
    SO_Cov(:,:)=x1;
    %% Kalman Gain
    K(1:9,1:6,t)=SO_Cov(:,:)*(InoCov(:,:))^-1;
    %% Mean State
    mu_t(:,t)=muBar_state(:,1,t)+K(:,:,t)*(real_obv(:,t)-muBar_obv(:,1,t));
    %% Covariance SP of the states
    covariance(:,:)=covariance(1:9,1:9)-K(:,:,t)*InoCov(:,:)*K(:,:,t)';
end

%% Plots
muBar_reducedobv(:,:)=muBar_obv(:,1,:);
t=1:1:numt;
figure(1), plot3(muBar_reducedobv(1,t),muBar_reducedobv(2,t),muBar_reducedobv(3,t),'.','DisplayName','Estimated Position');
hold on
plot3(real_obv(1,t),real_obv(2,t),real_obv(3,t),'DisplayName','Actual Position');
hold off;
title("Estimated Position vs Actual Trajectory")
xlabel('X component of Position');
ylabel('Y component of Position');
zlabel('Z component of Position');
legend('show');

figure(2), plot(1:numt,sqrt(sum(((muBar_reducedobv(1:3,t))-(real_obv(1:3,t))).^2)),'DisplayName','Estimated Position Error');
title("Position Error Convergence");
xlabel('Time (in Seconds)');
ylabel('Error between the measured and estimated values');
legend('show');

figure(3),plot3(muBar_reducedobv(4,t),muBar_reducedobv(5,t),muBar_reducedobv(6,t),'DisplayName','Estimated Magnetic Field');
hold on
plot3(real_obv(4,t),real_obv(5,t),real_obv(6,t),'r.','DisplayName','Actual Magnetic Field');
hold off;
title("Estimated Magnetic Field vs Actual Trajectory")
xlabel('X component of Magnetic field');
ylabel('Y component of Magnetic field');
zlabel('Z component of Magnetic field');
legend('show');
figure(4), plot(1:numt,sqrt(sum(((muBar_reducedobv(4:6,t))-(real_obv(4:6,t))).^2)),'--','DisplayName',' Estimated Magnetic Field Error');
title("Magnetic Field Error Convergence");
xlabel('Time (in Seconds)');
ylabel('Error between the measured and estimated values');
legend('show');



