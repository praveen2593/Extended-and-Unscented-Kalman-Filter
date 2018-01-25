%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function name: statevecInteg
% Returns state rates vector dX1 vector
% Input: 
% t = time
% s_vec=state vector,dX1 vector
% v= speed of the wire in m/s, scalar element
% omega= angular velocity of the magnetic field in rad/s, 3X1 vector
% P= position of the needle head (3X1 vector) from the s_vec
% H= Heading of the needle head (3X1 vector) from the s_vec
% B= Magnetic field applied to the needle head (3X1 vector) from the s_vec
% Output:
% s_vec_rate= state vector rate

% Name: Vineet Pandey
% CWID: 10826588
% Date: 11/22/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s_vec_rate] = statevecInteg(t,s_vec,v,omega)
P=s_vec(1:3);
H=s_vec(4:6);
B=s_vec(7:9);
e=1/3;
%% state rate estimation
s_vec_rate = zeros(9,1);
s_vec_rate(1:3)=H*v;
s_vec_rate(4:6)=-(cpMap(H)*cpMap(H))*B*v*e+H*(1-H'*H);
s_vec_rate(7:9)=-(cpMap(B))*omega+B*(1-B'*B);

end

