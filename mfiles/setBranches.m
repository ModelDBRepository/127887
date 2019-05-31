function [Np, L, Cm, Ra, Ena, Ek, El, r, gna, gk, gl, I0, xa, W, space_discr] = setBranches( numsections, Np_global, space_discr_global, inputflagpoint );


if numsections==1,

    k = 1;
    Np(k) = Np_global; % full size including ends
    % discretization method
    space_discr(k) = space_discr_global; % 0: cheby; 2,4,6: FD
    % constants
    L(k) = 400.*1e-4; % cm
    Cm(k) = 1; % microF/cm^2
    Ra(k) = .0354;
    Ena(k) = 50; % mV
    Ek(k) = -77;
    El(k) = -54.3;
    r(k) = 2*1e-4;
    % conductances set to passive
    gna(k) = 0; % 120 mS/cm^2
    gk(k) =  0;% 36; 
    gl(k) =  0.3;
    % input
    if 1==inputflagpoint,
        %NARROW
        I0(k) = (.65).*1e-3; % total current in microAmp
        middleofinput = ( 7 ./10).*L(k); % location
        W(k) = 1e-9.* (1 ./ 10).*L(k); % width
        xa(k) = middleofinput - W(k)./2;
    else
        % WIDE
        I0(k) = (.65).*1e-3; % total current in microAmp
        middleofinput = ( 5  ./10).*L(k); % location
        W(k) = ((9.999999999999)./10).*L(k); % width
        xa(k) = middleofinput - W(k)./2;
    end

else
    fprintf('*** ERROR: [%s] This file supports one section only\n',mfilename);
    return
end

