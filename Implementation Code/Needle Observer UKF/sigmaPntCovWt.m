%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function name: sigmaPntCovWt
% Returns two vectors with weights assigned to each sigma point and
% covariance around each point
% Input: 
% d = dimension of the state vector
% alpha=Factor deciding the spread of sigma points around the mean
% k= parameter that changes with the distribution of the data.k=2 for
% Gaussian distribution
% lambda = controls the spread of sigma points based on different distributions
% and dimensionality of state vectors
% Output:
% weights_pnt= weights corresponding to each sigma point,(2Xd+1,1)
% weights_cov= weights corresponding to covraiance,(2Xd+1,1)

% Name: Vineet Pandey
% CWID: 10826588
% Date: 11/22/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [weights_pnt,weights_cov] = sigmaPntCovWt(d)
alpha=0.7;
k=2;
beta=2;
weights_pnt=zeros(2*d+1,1);
weights_cov=zeros(2*d+1,1);
lambda=alpha^2*(d+k)-d;
for i=1:2*d+1
    if(i==1)
        weights_pnt(i)= lambda/(d+lambda);
    else
        weights_pnt(i)=1/(2*(d+lambda));
    end
end
for j=1:2*d+1
    if (j==1)
        weights_cov(j)=weights_pnt(1)+(1-alpha^2+beta);
    else
        weights_cov(j)=1/(2*(d+lambda));
        
    end
end
end

