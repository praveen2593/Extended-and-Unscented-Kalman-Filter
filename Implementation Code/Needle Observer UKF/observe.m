%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function name: observe
%Returns the output of the system given the states Position, P and Magnetic
%field, B
%[observation1] = observe(p,b)
%p = 3X3 position vector
%b = 3X3 magnetic field
%Cmat= Observation matrix:Identity matrix 6X6
%Observation1= Control Outputs 6X1 matrix
%Name: Vineet Pandey
%CWID: 10826588
%Date: 11/22/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [observation1] = observe(p,b)
observation1=zeros(6,1);
Cmat=[p;b];
observation1=eye(6)*Cmat;
end

