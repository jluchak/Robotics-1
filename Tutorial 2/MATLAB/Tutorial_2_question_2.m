%% Clear Workspace and Command Window
clc
clear

%% Given Constants and Variables
m=1.0; % Mass of link [kg]
l=0.2; % Length of link [m]
l_m=0.1; % Distance from the joint to the centre of mass of link [m]
I=0.01; % Moment of Inertia [kg*m^2]
T=2; % Total movement duration [s]
ts=0.01; % Timestep
t=0:ts:T; % time vector
w=t/T; % Dimensionless time
K=0.01; % Control Gain [s]
k=100; % Control Gain [Nm]

%% Given and determined Equations for desired position, velocity and acceleration of the Robot Endpoint
x_d=[0.273-0.2*(6*w.^5-15*w.^4+10*w.^3);0.273-0.1*(6*w.^5-15*w.^4+10*w.^3)]; % Position/trajectory of endpoint
x_d_dot=[-0.2/T^3*(30*t.^4/T^2-60*t.^3/T+30*t.^2);-0.1/T^3*(30*t.^4/T^2-60*t.^3/T+30*t.^2)]; % velocity of endpoint
x_d_dot_dot=[-0.2/T^3*(120*t.^3/T^2-180*t.^2/T+60*t.^1);-0.1/T^3*(120*t.^3/T^2-180*t.^2/T+60*t.^1)]; % Acceleration of endpoint

%% Determine Initial Angles of the desired and actual Shoulder Joints

x1=x_d(1); % x coordinate of endpoint
x2=x_d(2); % y coordinate of endpoint

q_d=zeros(2,T/ts+1); % Empty matrix for desired shoulder joint angles [left;right]
q_FB=zeros(2,T/ts+1); % Empty matrix for actual shoulder joints [left;right]
q_FF=zeros(2,T/ts+1); % Empty matrix for actual shoulder joints [left;right]

A=atan(x2(1)/x1(1)); % Intermediate step from geometry of question
B=acos(sqrt(x1(1)^2+x2(1)^2)/(2*l)); % Intermiediate step from geometry of question

q_d(1,1)=A+B; % Initial angle of left shoulder joint [rad]
q_d(2,1)=A-B; % Initial angle of right shoulder joint [rad]

q_FB(1,1)=q_d(1,1); % Initialize actual feedback left angle as equal to desired
q_FB(2,1)=q_d(2,1); % Initialize actual feedback right angle as equal to desired

q_FF(1,1)=q_d(1,1); % Initialize actual feedforward left angle as equal to desired
q_FF(2,1)=q_d(2,1); % Initialize actual feedforward right angle as equal to desired

%% Create and Initialize variables to find actual and desired q double dot (angular acceleration)

angular_acceleration_FB=zeros(2,T/ts+1); % Empty matrix for feedback angular acceleration
angular_velocity_FB=zeros(2,T/ts+1); % Empty matrix for feedback angular velocity
angular_acceleration_FF=zeros(2,T/ts+1); 
angular_velocity_FF=zeros(2,T/ts+1); 
angular_velocity_desired=zeros(2,T/ts+1);
angular_acceleration_desired=zeros(2,T/ts+1);

e=zeros(2,T/ts+1); % Empty matrix for Error function
e_dot=zeros(2,T/ts+1); % Empty matrix for derivative of error function

beta_FB=zeros(1,T/ts+1); % Empty vector for beta. Used in Matrix "H"
beta_FB(1)=2*m*l*l_m*cos(q_FB(2,1)-q_FB(1,1)); % initial value for beta
beta_FF=zeros(1,T/ts+1); % Empty vector for beta. Used in Matrix "H"
beta_FF(1)=2*m*l*l_m*cos(q_FF(2,1)-q_FF(1,1));
beta_desired=zeros(1,T/ts+1); % Empty vector for beta. Used in Matrix "H"
beta_desired(1)=2*m*l*l_m*cos(q_d(2,1)-q_d(1,1));

V_FB=zeros(2,T/ts+1); % Empty matrix for velocity vector
V_FF=zeros(2,T/ts+1);
V_desired=zeros(2,T/ts+1);

gamma_FB=zeros(1,T/ts+1); % Empty vector for gamma. Used in Matrix "V"
gamma_FB(1)=2*m*l*l_m*sin(q_FB(2,1)-q_FB(1,1)); % Initial value for gamma
gamma_FF=zeros(1,T/ts+1);
gamma_FF(1)=2*m*l*l_m*sin(q_FF(2,1)-q_FF(1,1));
gamma_desired=zeros(1,T/ts+1);
gamma_desired(1)=2*m*l*l_m*sin(q_d(2,1)-q_d(1,1));
%% Create a "for" loop to calculate angular position and trajectory of end point

for i=2:1:(T/ts+1)
    
    % Calculate and Print Desired angles
    
    q_d(1,i)=q_d(1,i-1)+angular_velocity_desired(1,i-1)*ts; % Desired Left shoulder angles printed into row 1 of matrix q_d
    q_d(2,i)=q_d(2,i-1)+angular_velocity_desired(2,i-1)*ts; % Desired right shoulder angles printed into row 2 of matrix q_d
    
    J=Jacobian2(q_d(1,i),q_d(2,i)); % Call the Jacobian FunctionJacobian2(q_d(1,2)
    angular_velocity_desired(:,i)=inv(J)*x_d_dot(:,i); % Desired angular velocity
    angular_acceleration_desired(:,i)=inv(J)*x_d_dot_dot(:,i); % Desired angular acceleration
    
    % Calculate and Print Actual angles
    q_FB(1,i)=q_FB(1,i-1)+angular_velocity_FB(1,i-1)*ts; % Actual left shoulder angles printed into row 1 of matrix q
    q_FB(2,i)=q_FB(2,i-1)+angular_velocity_FB(2,i-1)*ts; % Actual right shoulder angles printed into row 2 of matrix q
    
    e(:,i)=q_d(:,i)-q_FB(:,i); % Error function between desired and actual angles
    e_dot(:,i)=(e(:,i)-e(:,i-1))/ts; % Derivative of Error function
    
    torque_FB(:,i-1)=K*(e(:,i)+(k*e_dot(:,i))); % Linear Feedback controller torque
        
    alpha=2*m*l_m^2+m*l^2+2*I; % Intermediate equation "alpha" used in matrix H
    beta_FB(i)=2*m*l*l_m*cos(q_FB(2,i)-q_FB(1,i)); % Intermediate equation "beta" used in matrix H
    
    H=[alpha beta_FB(i);beta_FB(i) alpha]; % Mass matrix "H" of parallel robot
    
    gamma_FB(i)=2*m*l*l_m*sin(q_FB(2,i)-q_FB(1,i)); % Intermediate equation "gamma" used in Matrix V
    
    V_FB(:,i)=[0 -gamma_FB(i);gamma_FB(i) 0]*[angular_velocity_FB(1,i)^2;angular_velocity_FB(2,i)^2]; % Velocity Matrix
    
    angular_acceleration_FB(:,i)=inv(H)*(torque_FB(:,i-1)-V_FB(:,i)); % Calculate angular acceleration
    angular_velocity_FB(:,i)=angular_velocity_FB(:,i-1)+angular_acceleration_FB(:,i-1)*ts; % Calculate angular velocity
end  
%% Feedforward Controller
for i=2:1:(T/ts+1)
    
    q_FF(1,i)=q_FF(1,i-1)+angular_velocity_FF(1,i-1)*ts; % Actual left shoulder angles printed into row 1 of matrix q
    q_FF(2,i)=q_FF(2,i-1)+angular_velocity_FF(2,i-1)*ts; % Actual right shoulder angles printed into row 2 of matrix q
     
    % Useing desired angles to find feedforward torque   
    beta_desired(i)=2*m*l*l_m*cos(q_d(2,i)-q_d(1,i)); % Intermediate equation "beta" used in matrix H
    H_desired=[alpha beta_desired(i);beta_desired(i) alpha]; % Mass matrix "H" of parallel robot
    gamma_desired(i)=2*m*l*l_m*sin(q_d(2,i)-q_d(1,i)); % Intermediate equation "gamma" used in Matrix V
    V_desired(:,i)=[0 -gamma_desired(i);gamma_desired(i) 0]*[angular_velocity_desired(1,i)^2;angular_velocity_desired(2,i)^2]; 
    % Velocity Matrix
     
    % Calculate Feedforward torque and superpositioned torque
    torque_FF(:,i-1)=H_desired*angular_acceleration_desired(:,i)+V_desired(:,i);
    torque(:,i-1)=torque_FF(:,i-1)+torque_FB(:,i-1); 
    
    % Calculate feedforward angles
    beta_FF(i)=2*m*l*l_m*cos(q_FF(2,i)-q_FF(1,i)); % Intermediate equation "beta" used in matrix H
    H_FF=[alpha beta_FF(i);beta_FF(i) alpha]; % Mass matrix "H" of parallel robot
    gamma_FF(i)=2*m*l*l_m*sin(q_FF(2,i)-q_FF(1,i)); % Intermediate equation "gamma" used in Matrix V
    V_FF(:,i)=[0 -gamma_FF(i);gamma_FF(i) 0]*[angular_velocity_FF(1,i)^2;angular_velocity_FF(2,i)^2]; % Velocity Matrix
           
    angular_acceleration_FF(:,i)=inv(H_FF)*(torque(:,i-1)-V_FF(:,i)); % Calculate angular acceleration
    angular_velocity_FF(:,i)=angular_velocity_FF(:,i-1)+angular_acceleration_FF(:,i-1)*ts; % Calculate angular velocity
end
%% Calculate Trajectory of Endpoint
% Feedback actual
FB_actual_x_traj=l*cos(q_FB(2,:))+l*cos(q_FB(1,:));
FB_actual_y_traj=l*sin(q_FB(2,:))+l*sin(q_FB(1,:));
% Feedforward actual
FF_actual_x_traj=l*cos(q_FF(2,:))+l*cos(q_FF(1,:));
FF_actual_y_traj=l*sin(q_FF(2,:))+l*sin(q_FF(1,:));

%% Change calculated angle from radians to degrees
% Desired 
Desired_L =q_d(1,:)*180/pi();
Desired_R =q_d(2,:)*180/pi();
% Feedback Actual
Actual__FB_L =q_FB(1,:)*180/pi();
Actual__FB_R =q_FB(2,:)*180/pi();
% Feedforward Actual
Actual_FF_L =q_FF(1,:)*180/pi();
Actual_FF_R =q_FF(2,:)*180/pi();

%% Plot Desired and Actual Angles against time
figure(1)
% ql
subplot(2,2,1:2)
hold on
plot(t,Desired_L,'b'); % Desired plot
plot(t,Actual__FB_L,'r');
plot(t,Actual_FF_L,'g--');
legend('Desired','Actual FB','Actual FF','location','northwest');
title('Left Shoulder "ql"') % Subplot title
xlabel('Time [s]');
ylabel('Angle [deg]');
hold off

% qr
subplot(2,2,3:4)
hold on
plot(t,Desired_R,'b'); % Desired right plot
plot(t,Actual__FB_R,'r'); % Actual right plot
plot(t,Actual_FF_R,'g--');
legend('Desired','Actual FB','Actual FF','location','northeast');
title('Right Shoulder "qr"') % Subplot title
xlabel('Time [s]');
ylabel('Angle [deg]');
hold off

%% Plot Desired and Actual endpoint positions
figure(2)
% X position against time
subplot(2,2,1:2)
hold on
plot(t,x_d(1,:),'b'); % Desired
plot(t,FB_actual_x_traj,'r'); % Actual feedback
plot(t,FF_actual_x_traj,'g--'); % Actual feedforward
legend('Desired','Actual FB','Actual FF','location','northeast');
xlabel('Time [s]'); % Time in seconds
ylabel('Position "x1" [m]'); % x coordinate
hold off

% Y position against time
subplot(2,2,3:4)
hold on
plot(t,x_d(2,:),'b'); % Desired
plot(t,FB_actual_y_traj,'r'); % Actual feedback
plot(t,FF_actual_y_traj,'g--'); % Actual feedforward
legend('Desired','Actual FB','Actual FF','location','northeast');
xlabel('Time [s]'); % Time in seconds
ylabel('Position "x2" [m]'); % y coordinate
hold off

%% Plot Desired and Actual Endpoint Trajectories
figure(3)
hold on
plot(x_d(1,:),x_d(2,:),'b'); % Desired
plot(FB_actual_x_traj,FB_actual_y_traj,'r'); % Actual feedback
plot(FF_actual_x_traj,FF_actual_y_traj,'g--'); % Actual feedforward
legend('Desired','Actual FB','Actual FF','location','southeast');
title('Desired and Actual Endpoint Trajectory of Robot');
xlabel('Position "x1" [m]'); % x coordinate
ylabel('Position "x2" [m]'); % y coordinate
hold off