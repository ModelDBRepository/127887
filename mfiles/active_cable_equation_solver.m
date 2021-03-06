%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Active branching cable neuronal simulator.
%  Compares numerical with exact solution and plots the results for
%  multiple simulations.
%  Space discretization methods: Chebyshev spectral or Finite Difference
%  (FD) 2nd-6th orders.
%  (C) Ahmet Omurtag 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET PARAMETER VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tree params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of sections
numsections = 1;

% Tree structure as boundary conditions
[BC_branchends, BC_nodecurrents, BC_nodeequality ] = setNodes(numsections);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accuracy of time stepping
ATOL = 3e-12;
RTOL = ATOL;
fprintf('Accuracy of time integration is set to Relative %e, Absolute %e\n',RTOL,ATOL);

% Total simulation time
T = 10; % ms

% Input current start and end times
t_input = [0, 1000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop params for multiple simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop twice (point input and broad input)
for  inputflagpoint = [0 1],

    % Space discretization methods; 
    % space_discr_global = 
    % 0: Spectral
    % 2: Second order Finite Difference
    % 4: Fourth order Finite Difference
    % 6: Sixth order Finite Difference

    countS = 1;
    for space_discr_global = [0 2 4 6],

        % Grid sizes Np (includes end points)
        countN = 1;
        for Np_global = [ 8 12 16 24 32 48 64 ]

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set the individual section params
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Set parameters and synaptic input for each section
            [Np, L, Cm, Ra, Ena, Ek, El, r, gna, gk, gl, I0, xa, W, space_discr] = setBranches( numsections, Np_global, space_discr_global, inputflagpoint );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SET UP CALCULATIONS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Preliminary calculations for each section
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Sections
            for k = 1:numsections,
                [F,G,x] = getDifferentiationMatrix( space_discr(k), Np(k), L(k) );

                % Create discretized operators
                P = getLHSOperator( G, Np(k), 0 );
                Q0 = getStructureMatrix( P, Np(k), r(k).*ones(Np(k),1), gl(k).*ones(Np(k),1), Ra(k), F );

                % Right hand side vectors related to synaptic input
                R0 = getRHS( P, El(k), gl(k).*ones(Np(k),1) );
                J = getRHSInput( I0(k), r(k), W(k), xa(k), Np(k), x, 0 );
                P = P.*Cm(k); % factor in capacitance

                label = ['sect' num2str(k)];
                if k==1,
                    FS = struct(label,F);
                    PS = struct(label,P);
                    Q0S = struct(label,Q0);
                    R0S = struct(label,R0);
                    JS = struct(label,J);
                    xS = struct(label,x);
                else
                    FS = setfield(FS,label,F);
                    PS = setfield(PS,label,P);
                    Q0S = setfield(Q0S,label,Q0);
                    R0S = setfield(R0S,label,R0);
                    JS = setfield(JS,label,J);
                    xS = setfield(xS,label,x);
                end
            end
            clear F G P Q0 R0 J

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Tree calculations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Boundary condition matrix
            [BC,ColIndx] = calcBCMatrix(numsections,Np,r(1,:),r(end,:),Ra,FS,BC_branchends,BC_nodecurrents,BC_nodeequality);
            K = getK(BC,ColIndx); % K converts v to u, i.e. u = K*v

            % Tree operators (stack section operators diagonally into the big operator)
            PP = calcTreeOper(numsections,Np,PS);
            QQ0 = calcTreeOper(numsections,Np,Q0S);
            RR0 = calcTreeRHS(numsections,Np,R0S);
            JJ = calcTreeRHS(numsections,Np,JS);

            % Reduce LHS by using BC matrix
            PPbar = reduceByBC(PP,K,ColIndx);

            % Generate total LHS oper (for system including channel vars)
            PPvmhn = calcBigLHSOper(PPbar,Np);

            % Number of membrane potential variables
            Nv = sum(Np) - numsections.*2;

            % Num of each type of channel var
            Nm = sum(Np);

            % Full conductance utilities for vector field
            ggna = []; ggk = [];
            for k = 1:numsections,
                ggna = [ ggna; gna(k).*ones(Np(k),1) ];
                ggk = [ ggk; gk(k).*ones(Np(k),1) ];
            end
            ggEna = []; ggEk = [];
            for k = 1:numsections,
                ggEna = [ ggEna; Ena(k).*gna(k).*ones(Np(k),1) ];
                ggEk = [ ggEk; Ek(k).*gk(k).*ones(Np(k),1) ];
            end
            clear PS Q0S R0S JS

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Initialize dependent variables 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            vini = El(1); % Set initial membrane potential to leakage reversal potential
            mhnini = 0;
            vv = vini.*ones(Nv,1);
            mhn = mhnini.*ones(3.*Nm,1); % set channel variables to unity
            vmhn = [ vv; mhn ];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TIME STEPPING
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            options = odeset('RelTol',RTOL,'AbsTol',ATOL,'Mass',PPvmhn);
            Tspan = [0 T];

            [TT,YY] = ode15s( @(t,y) cable_vectfield_active5(t,y,numsections,t_input,Np,PP,QQ0,RR0,JJ,K,ColIndx,ggna,ggk,ggEna,ggEk,PPbar) ...
                ,Tspan,vmhn,options );
            
            % Store result
            VV = YY(:,1:Nv);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Membrane potential: numerical and exact solutions, passive, one section
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if numsections == 1,
                figure(1); clf
                xx = x(2:end-1);
                % Downsample for plotting
                if length(TT) > 50,
                    numtimesteps = fix(length(TT)./50); TTn = TT(1:numtimesteps:end); VVn = VV(1:numtimesteps:end,:);
                else
                    TTn = TT; VVn = VV;
                end
                mesh(xx,TTn,VVn); xlabel('x (cm)'); ylabel('Time (ms)'); zlabel('V (mV)');
                title('Space-time plot of the membrane potential');
            else
                fprintf('*** ERROR: [%s] This file supports one section only\n',mfilename);
                return
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate error at each grid point
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            xx = x(2:end-1);
            Vexact = cable_exact3(xa,W,r,L,Ra,0,gl,xx,I0,Cm,T) + El;
            Err = Vexact - VV(end,:)';

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Final membrane potential and exact solution for passive one section
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if numsections == 1,
                figure(2);clf;
                % Need exact soln for plot
                deltt = L./50;
                xxplot = 0:deltt:L;
                Vexactplot = cable_exact3(xa,W,r,L,Ra,0,gl,xxplot',I0,Cm,T) + El;
                plot(xxplot,Vexactplot,xx,VV(end,:),'r.'); xlabel('x (cm)'); ylabel('V(x,T)');
                legend('Exact solution','Numerical solution',0); legend boxoff; title('Comparison of numerical and exact solutions');
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % End grid size loop
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            err =  mean(abs(Err));;
            fprintf('Space discr method: %d, Grid size  N=%d, Mean Err %g\n',space_discr_global,Np_global, err);
            na = max(find(x<xa));
            ErrorV(countS,countN) = err;
            NN(countN) = Np(1);
            countN = countN + 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % End Space discretization loop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        countS = countS + 1;
    end

    if 1==inputflagpoint, ErrorVPt = ErrorV; else, ErrorVDist = ErrorV; end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACCURACY RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot numerical errors as function of grid size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numsections==1,
    figure(100); clf
    % Main plot
    % Point input result
    loglog(NN,ErrorVPt(1,:),'r.--'); hold on;
    loglog(NN,ErrorVPt(2,:),'b.--'); hold on;
    loglog(NN,ErrorVPt(3,:),'g.--'); hold on;
    loglog(NN,ErrorVPt(4,:),'k.--'); hold on;
    % Distributed input result
    p1=loglog(NN,ErrorVDist(1,:),'r.-'); hold on;
    p2=loglog(NN,ErrorVDist(2,:),'b.-'); hold on;
    p3=loglog(NN,ErrorVDist(3,:),'g.-'); hold on;
    p4=loglog(NN,ErrorVDist(4,:),'k.-'); hold on;    
    xlabel('N','fontsize',22)
    ylabel('E_N','fontsize',22)
    axis([NN(1)  NN(end) min(min([ErrorVDist; ErrorVPt]))  1]);
    legend([p1 p2 p3 p4],'Spectral','FD2','FD4','FD6','location','southwest'); legend boxoff
end


