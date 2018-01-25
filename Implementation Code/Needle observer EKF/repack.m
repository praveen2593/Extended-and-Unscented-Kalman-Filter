function [X] = repack(P,Hs,B,Jx,Ju)
if ~isempty(Jx)
    X = zeros(126,1);
else
    X = zeros(9,1);
end

X(1:3) = P;
X(4:6) = Hs;
X(7:9) = B;
if ~isempty(Jx)
X(10:90) = reshape(Jx,[81,1]);
X(91:end) =reshape(Ju,[36,1]);
end
end