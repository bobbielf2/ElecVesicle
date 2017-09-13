function Main(RunFile)

savedata = 0; % save data?

% Preliminary
% -------------------------------------------------------------------------
warning off
if nargin==0
    run('Submit')
else
    run(RunFile)
end
X=cell(1,2);
Y=X;
SIG=cell(1,1);
M=length(y);
LS=zeros(M,1);
for k=1:M
    LS(k)=length(x{k})-1;
end
X{1}=x;
Y{1}=y;
G=[];
G1=[];
ArcL=zeros(M,1);
for k=1:M
    ArcL(k)=ArcLength(x{k},y{k},2*pi);
end
AreaA=zeros(M,1);
for k=1:M
    AreaA(k)=Area(x{k},y{k});
end
if ~exist('Name','var')
    Name='Test'; % Name of output file
end
Time=zeros(T,1);

%%%%%%%%%%%% EHD Param %%%%%%%%
f = cell(M,1);
F_el = f;
Vm = zeros(LS(1),M); % jump of electric potential
sigma_i = conduct_rat*ones(M,1); % conductivity interior
sigma_e = 1; % conductivity exterior
Gm = Gm*ones(M,1); % conductance
Cm = 1*ones(M,1); % capacitance
%%%%%%%%%%%% EHD Param %%%%%%%%
% -------------------------------------------------------------------------


% Construct arc length correction operator
ALC=@(alpha,beta)ALC_operator(alpha,beta,ArcL,dt);

% Reparameterize and initial save
[X,Y,SIG]=Reparameterize(X,Y,SIG,1);
X_1=X{1};
Y_1=Y{1};
if savedata
  save(Name,'X_1','Y_1','dt','kb','mu','k','E_infty')
end

ss = cell(M,1);
for j = 1:M
    ss{j}.x = X{1}{j}(2:end)+1i*Y{1}{j}(2:end);
    ss{j} = quadr(ss{j});
end
[SLP,SLPprime] = SLP_VV_lap(ss);
L_Op = L_operator_Vm(SLPprime,sigma_i,sigma_e);
    
    
tic
for k=1:T
    
    % Initial guess for GMRES
    if k>1
        if k>2
            G1=2*G-G0;
        else
            G1=G;
        end
    end
    G0=G;
    
    n=min(k,2);
    
    % Reparameterize x_k, y_k, and sig_k with respect to arc length
    [X,Y,SIG]=Reparameterize(X,Y,SIG,n);
    
    % Apply filter to x_k and y_k
    [X,Y]=Filter(X,Y,n);
    
    % Construct Fb operator
    Fb=@(alpha,beta)Fb_operator(alpha,beta,X{n},Y{n},kb);
    
    % Construct Fsig operator
    Fsig=@(gamma)Fsig_operator(gamma,X{n},Y{n});
    
    % Construct Vesicle-Vesicle interaction operator
        % Reformating
        s=cell(M,1);
        for j=1:M
            s{j}=[X{n}{j}(2:end),Y{n}{j}(2:end)];
        end
        sf=zeros(2*sum(LS),1);
        ls=0;
        for j=1:M
            sf(ls+1:ls+LS(j))=X{n}{j}(2:end);
            sf(end/2+(ls+1:ls+LS(j)))=Y{n}{j}(2:end);
            ls=ls+LS(j);
        end
    
        % Self interactions
        SS=SLP_self_operator(X{n},Y{n});
        
        % Close interactions
        CM=IntMat(sf,LS,sf);
    
        % Far interactions
        FM=MakeFM(s);
    
        % Create interaction operator
        SLP_VV=@(f)SLP_VV_operator2(f,X{n},Y{n},SS,CM,FM,mu);
    
    % compute electrical force
    ss = cell(M,1);
    for j = 1:M
        ss{j}.x = X{n}{j}(2:end)+1i*Y{n}{j}(2:end);
        ss{j} = quadr(ss{j});
    end
    q = get_q(ss,Vm,SLP,L_Op,E_infty,sigma_i,sigma_e,tol);
    tau = f_el(ss, SLP, SLPprime, E_infty, Vm, q);
    for j=1:M
        F_el{j}=[tau([end/2,1:end/2],j),tau([end,end/2+1:end],j)];
    end
        
    % Construct inextensibility operator
    P=@(alpha,beta)P_operator([alpha;beta],X{n},Y{n});
    
    % construct v_inf operator
    v_inf=@(alpha,beta)v_inf_operator(alpha,beta,X{n},shear);
    
    % Construct Stokes Linear Operator
    Fdt=@(alpha,beta,gamma)Fdt_operator(alpha,beta,gamma,Fb,Fsig,dt);
    %LinOP=@(x)[P(x(1:end/3),x(end/3+1:2*end/3));x(1:2*end/3)-SLP_VV(Fdt(x(1:end/3),x(end/3+1:2*end/3),x(2*end/3+1:end)))];
    LinOP=@(x)[P(x(1:end/3),x(end/3+1:2*end/3));x(1:2*end/3)-SLP_VV(Fdt(x(1:end/3),x(end/3+1:2*end/3),x(2*end/3+1:end)))-dt*v_inf(x(1:end/3),x(end/3+1:2*end/3))];
    
    % Create Preconditioner
    PreInv=PreInv_operator(X{n},Y{n},SS,kb,dt,mu);
    PreCond=@(x)PreCond_operator(x(1:end/3),x(end/3+1:2*end/3),x(2*end/3+1:end),PreInv,M,LS);
    
    % Compute bending force
    fb=Fb(X{n},Y{n});
    
    % Bending & electric (for constructing RHS)
    for i = 1:M
        f{i} = fb{i}-F_el{i};
    end
    
    % Compute RHS vector
    %rhs=[ALC(X{n},Y{n});SLP_VV(fb)];
    rhs=[ALC(X{n},Y{n});SLP_VV(f)+V_inf(X{n},Y{n},shear)];
    
    % Solve for u_x, u_y, and sig using GMRES
    G=gmres(LinOP,rhs,[],tol,20,PreCond,[],G1);
    
    X{1}=X{n};
    Y{1}=Y{n};
    
    % Update results
    ls=0;
    for j=1:M
        X{2}{j}=X{1}{j}+dt*G([ls+1:ls+LS(j),ls+1]);
        Y{2}{j}=Y{1}{j}+dt*G([end/3+ls+1:end/3+ls+LS(j),end/3+ls+1]);
        SIG{1}{j}=G([2*end/3+ls+1:2*end/3+ls+LS(j),2*end/3+ls+1]);
        ls=ls+LS(j);
    end
    
    
    % update membrane potential Vm
      % construct lin op for evolving Vm
      [SLP,SLPprime] = SLP_VV_lap(ss);
      L_Op = L_operator_Vm(SLPprime,sigma_i,sigma_e);
      LinOPVm = @(alpha) LinOP_Vm(alpha,ss,SLP,L_Op,Cm,Gm,dt,sigma_i, sigma_e);
      
      % construct rhs
      N = length(ss{1}.x);
      CM = kron(diag(Cm),eye(N)); % Cm matrix
      lambda = sigma_i*sigma_e./(sigma_i+sigma_e);
      Lambda = kron(diag(lambda),eye(N)); % Lambda matrix
      Nx = []; % outward normals on bdry
      for i = 1:M, Nx = [Nx;ss{i}.nx]; end
      rhs_elec = CM*L_Op*Vm(:) + dt*Lambda*real(conj(E_infty)*Nx);
      G2 = gmres(LinOPVm,rhs_elec,[],tol,length(rhs_elec));
      Vm = reshape(G2,N,M);
    
    
    % Apply area correction
    for j=1:M
        c1=-Area(Dtp(Y{2}{j}),Dtp(X{2}{j}));
        c2=Area(Dtp(Y{2}{j}),Y{2}{j})-Area(X{2}{j},Dtp(X{2}{j}));
        c3=Area(X{2}{j},Y{2}{j})-AreaA(j);
        disc=sqrt(c2^2-4*c1*c3);
        if abs((-c2+disc)/(2*c1))>abs((-c2-disc)/(2*c1))
            a=(-c2-disc)/(2*c1);
        else
            a=(-c2+disc)/(2*c1);
        end
        X{2}{j}=X{2}{j}+a*Dtp(Y{2}{j});
        Y{2}{j}=Y{2}{j}-a*Dtp(X{2}{j});
    end
    
    % Save Results
    if savedata
      eval(['X_',num2str(k+1),'=X{2};']);
      eval(['Y_',num2str(k+1),'=Y{2};']);
      eval(['SIG_',num2str(k),'=SIG{1};']);
      %%%%%%%%%%%%%
      eval(['F_el_',num2str(k),'=F_el;']);
      eval(['F_',num2str(k),'=f;']);
      eval(['Vm_',num2str(k),'=Vm;']);
      eval(['q_',num2str(k),'=q;']);
      eval(['phi_in_',num2str(k),'=phi_in;']);
      eval(['phi_ex_',num2str(k),'=phi_ex;']);
      %%%%%%%%%%%%%
      save(Name,['X_',num2str(k+1)],['Y_',num2str(k+1)],...
          ['SIG_',num2str(k)],['F_el_',num2str(k)],['F_',num2str(k)],...
          ['Vm_',num2str(k)], ['q_',num2str(k)],...
          ['phi_in_',num2str(k)],['phi_ex_',num2str(k)],'k','G','G0','-append')
      % free up memory
      eval(['clearvars X_',num2str(k+1), ' Y_',num2str(k+1), ' SIG_',num2str(k), ' F_el_',num2str(k)])
    end
    
    % Estimates time remaining
    Time(k)=toc;tic
    Tav=sum(Time(1:k))/k;
    ET=(T-k)*Tav;
    disp([num2str(round(100*k/T)),...
        '% Complete / Estimated time remaining: ',...
        num2str(floor(ET/3600)),' hours, ',...
        num2str(round(floor(ET/60)-60*floor(ET/3600))),' minutes, ',...
        num2str(round(floor(ET)-60*floor(ET/60))),' seconds'])
    
    if k>1
        vel_max = 0;
        figure(1);
        for i = 1:M
            plot(X{2}{i},Y{2}{i},'r')
            hold on
            %plot(X{2}{i}(1),Y{2}{i}(1),'ro','Markersize',10) %track rotation
            Vel = [X{2}{i}-X{1}{i}, Y{2}{i}-Y{1}{i}];
            vel_temp = max(max(abs(Vel)));
            %quiver(X{2}{i},Y{2}{i}, Vel(:,1), Vel(:,2)); %plot velocity
            if vel_max<vel_temp, vel_max = vel_temp; end
        end
        disp('max vel: ')
        disp(vel_temp)
        axis([-2,2, -2,2]/.6)
        axis equal
        title(num2str(k))
        %pbaspect([2,1,1])
        drawnow
        %F(k) = getframe(h);
        hold off
        if vel_max>10
            disp('Vel too big!!')
            break
        end
    end
end
