function xDot = needleModel(t,x,u,c)
% x =  [Shaft Pos; Shaft Heading; B-Field Direction; Unpacked Jx; unpacked
% Ju]
% u = [insertion velocity, B anglar Velocity
% t = time
% c = Curvature Constant relating missalignment between head and shaft with
%    resulting path curvature
if ~exist('c','var')
c = 7*sqrt(2)/2 * 1/3;
end
[p,hs,B,Jx,Ju] = unpack(x);

v = u(1);
w = u(2:4);
I = eye(3);

pDot = hs*v;
hsDot = v*c*(I-hs*hs')*B+hs*(1-hs'*hs);
BDot = cross(w,B)+B*(1-B'*B);
if length(x) > 9
    df_dx = [zeros(3,3) I*v zeros(3,3);
        zeros(3,3) -c*v*(hs'*B*eye(3)+hs*B')+(1-hs'*hs)*eye(3)-2*(hs*hs') c*v*(I-hs*hs');
        zeros(3,6) cpMap(w)+I*(1-B'*B)-2*(B*B')];
    
    Jx_dot = df_dx*Jx;
    
    Ju_dot = df_dx*Ju + [[hs;c*(I-hs*hs')*B;zeros(3,1)] [zeros(6,3);-cpMap(B)]];
else
    Jx_dot = [];
    Ju_dot = [];
end
xDot = repack(pDot,hsDot,BDot,Jx_dot,Ju_dot);

end

function xCross = cpMap(x)
xCross = [0 -x(3) x(2);
    x(3), 0 , -x(1);
    -x(2) x(1) 0];
end

