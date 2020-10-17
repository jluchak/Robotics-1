%% Constant Variables
T=2; % Total time duration [s]
ts=0.01; % time step
d=0; % Distance between the shoulders
a=0.3; % Length of upper arm [m]
b=0.3; % Length of forearm [m]
c=0.2; % Length of hand [m]
ql=deg2rad(150); % Initial left motor angle
qr=deg2rad(30); % Initial right motor angle
theta=deg2rad(45); % Initial hand orientation
%% Initiate angular joint matrix
joint_angles=zeros(3,T/ts); % Angle profiles
joint_angles(1,1)=ql;
joint_angles(2,1)=qr;
joint_angles(3,1)=theta;
%% Define initial and final position of point H
x0=[0.141,0.441]'; % Initial position
x1=[0.241,0.641]'; % Final position
length=sqrt((x1(2)-x0(2))^2+(x1(1)-x0(1))^2); % length of vector
%% Define Velocity profile
t=0:ts:T; % time matrix
tau=t/T; % non-dimensional time
sigma=30*tau.^2.*(tau.^2-2.*tau+1); % Velocity Profile
Real_Velocity=sigma*length; % Real Velocity
%% Define Displacement/Position profiles
posx=0:1:T/ts; % Create position x matrix
posy=0:1:T/ts; % Create position y matrix
posx(1,1)=x0(1); % Set starting position x
posy(1,1)=x0(2); % Set starting position y
sum=0:1:T/ts; % Create matrix for displacement
%% Calculate Position and Trajectory
for i=1:1:T/ts
    sum(1,i+1)=sum(1,i)+ts*(length)*sigma(1,i)/T; % calculate displacements
    posx(1,i+1)=posx(1,i)+ts*(length)*sigma(1,i)/T*(x1(1)-x0(1))/(length); % calculate position x
    posy(1,i+1)=posy(1,i)+ts*(length)*sigma(1,i)/T*(x1(2)-x0(2))/(length); % calculate position y
end
%% Calculate Joint Angles: ql, qr, theta 
joint_velocity_norm=zeros(1,T/ts+1); % Set up an emprty matrix to store joint veloicty norm values
for i=1:1:T/ts    
J=Jacobian(joint_angles(1,i),joint_angles(2,i),joint_angles(3,i)); % Call the Jacobian Function
J_inv=J.'*inv(J*J.'); % Calculate pseudo inverse Jacobian
velocity=J_inv*[sigma(1,i+1)*length*(x1(1)-x0(1))/length;sigma(1,i+1)*length*(x1(2)-x0(2))/length];

joint_angles(1,i+1)=joint_angles(1,i)+ts/T*velocity(1,1); % ql
joint_angles(2,i+1)=joint_angles(2,i)+ts/T*velocity(2,1); % qr
joint_angles(3,i+1)=joint_angles(3,i)+ts/T*velocity(3,1); % theta

joint_velocity_norm(i)=norm(velocity); % Store the norm values
end
%% Change joint angles from radians to degrees
Deg_joint_angles=joint_angles*180/pi();
%% Print Plots for Question 2B
figure(1)
hold on
% Hand velocity
subplot(3,2,1:2)
plot(tau*T,Real_Velocity,'b-');
title('Hand Velocity Profile');
xlabel('Time [s]');
ylabel('Velocity [m/s]');
% Hand trajectory
subplot(3,2,3:4)
plot(posx,posy);
title('Hand Trajectory Profile')
xlabel('Position "x1" [m]');
ylabel('Position "x2" [m]');
axis([0.13 0.25 0.4 0.7]);
% Position Profile
subplot(3,2,5:6)
plot(t,sum);
title('Hand Position Profile');
xlabel('Time [s]');
ylabel('Displacement [m]');
hold off
%% Print plots for question 2C
% ql versus time
figure(2)
hold on
subplot(3,2,1:2)
plot(t, Deg_joint_angles(1,:));
title('Left Motor Angle "ql(t)" Profile');
xlabel('Time [s]');
ylabel('Joint Angle [Deg]');
% qr versus time
subplot(3,2,3:4)
plot(t, Deg_joint_angles(2,:));
title('Right Motor Angle "qr(t)" Profile');
xlabel('Time [s]');
ylabel('Joint Angle [Deg]');
% theta versus time
subplot(3,2,5:6)
plot(t, Deg_joint_angles(3,:));
title('Hand Orientation Angle "theta(t)" Profile');
xlabel('Time [s]');
ylabel('Joint Angle [Deg]');
hold off
% Plot the joint velocity norm
figure(3)
plot(t,joint_velocity_norm);
title('Joint Velocity Norm Versus Time');
xlabel('Time [s]');
ylabel('Joint Velocity Norm [deg/s]');

%% End of Script

