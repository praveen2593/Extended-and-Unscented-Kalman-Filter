function [P,Hs,B,Jx,Ju] = unpack(x)
P = reshape(x(1:3),[3,1]);
Hs = reshape(x(4:6),[3,1]);
B = reshape(x(7:9),[3,1]);

if length(x)>9
Jx = reshape(x(10:90),[9,9]);
Ju = reshape(x(91:end),[9,4]);
else
    Jx = [];
    Ju = [];
end

end

