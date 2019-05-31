function dvmhndt = cable_vectfield_active5( tt, yy ...
    , numsections, tt_input, Np, PP, QQ0, RR0, JJ, K, ColIndx, ggna, ggk, ggEna, ggEk, PPbar ...
    )

clear PPbar

Ntot = size(yy,1);
Nm = sum(Np);
Nv = Nm - size(ColIndx,2);

% extract mem pop and channel vars for full system
vv = yy(1:Nv); yy(1:Nv) = []; 
mm = yy(1:Nm); yy(1:Nm) = []; 
hh = yy(1:Nm); yy(1:Nm) = []; 
nn = yy(1:Nm); yy(1:Nm) = [];
m3h = mm.*mm.*mm.*hh;
n4 = nn.*nn.*nn.*nn;

% calc v operators
QQ = QQ0 - PP * diag( ggna.*m3h + ggk.*n4 );
QQbar = reduceByBC(QQ,K,ColIndx);
RR = RR0 + PP*( ggEna.*m3h + ggEk.*n4 );
if tt_input(1) <= tt & tt <= tt_input(2),
    RR = RR + JJ;
end
    
% vec field of mem pot
dvdt = QQbar*vv + RR;

% calc full mem pot (including ends)
uu = K*vv;
vvfull(ColIndx) = uu;
numcols = Nm;
NN = [1:numcols];
indx = setdiff(NN,ColIndx); % the reduced system indices
vvfull(indx) = vv;

% calc vec field of channel vars
[am,bm,ah,bh,an,bn] = rates(vvfull');

dmdt = am.*(1-mm) - bm.*mm;
dhdt = ah.*(1-hh) - bh.*hh;
dndt = an.*(1-nn) - bn.*nn;

clear am ah an bm bh bn

% combined vector field
dvmhndt = [
    dvdt; ...
    dmdt; ...
    dhdt; ...
    dndt ...
    ];
