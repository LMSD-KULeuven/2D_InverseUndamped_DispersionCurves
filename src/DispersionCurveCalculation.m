%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%        INVERSE APPROACH FOR DISPERSION CURVE CALCULATIONS        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT DATA (definition of the problem)
% UC size
Lx = 0.05;                           % UC size x-direction [m]
Ly = 0.05;                           % UC size y-direction [m]
Lz = 0.005;                          % UC size z-direction [m]
% Material
E = 210e9;                           % Young's modulus [Pa]
v = 0.3;                             % Poisson ratio [-]
rho = 7800;                          % Mass density [kg/m^3]
% Mesh definition
n_elx = 10;                          % Number of elements in x-dir [-] (>0)
n_ely = 10;                          % Number of elements in y-dir [-] (>0)
n_elz = 3;                           % Number of elements in z-dir [-] (>0)
% Addition to UC bare plate
scatterer = 'none';                  % 'none', 'resonator', 'mass'
f_res = 2500;                        % resonance of resonator [Hz]
m_ratio = 0.3;                       % Scatterer mass ratio vs plate [-]
% Parameters for dispersion curves
cont_name = {'O','A','B','O'};       % Name of the IBC contour
cont_co =  [0,0; pi,0; pi,pi; 0,0];  % IBC definition
step_size = 0.01*pi;                 % Resolution propagation constant
n_curves = 10;                       % Number of calculated wave modes 

%% FE preprocessing
% Element size
ax = Lx/n_elx;
ay = Ly/n_ely;
az = Lz/n_elz;
% Information UC
n_elem = n_elx*n_ely*n_elz;   
n_nodes = (n_elx+1)*(n_ely+1)*(n_elz+1);
nDOF = 3;                           % Number of DOFs per node;
n_DOFs = nDOF*n_nodes;
mass_UC = rho*(n_elx*ax)*(n_ely*ay)*(n_elz*az);
% Element matrices
[KE,ME] = KM(E,v,rho,ax,ay,az);
% FEM information (hard-coded for nDOF = 3)
nodeNrs = reshape(1:n_nodes, 1+n_ely, 1+n_elz, 1+n_elx); % nodes numbering
cMat = reshape(nDOF * nodeNrs(1:n_ely, 1:n_elz, 1:n_elx)+1, n_elem, 1) + ...
    [0,1,2,3*(n_ely+1)*(n_elz+1)+[0,1,2,-3,-2,-1],-3,-2,-1,3*(n_ely + ...
    1)+[0,1,2],3*(n_ely+1)*(n_elz+2)+[0,1,2,-3,-2,-1],3*(n_ely+1)+[-3,-2,-1]]; % connectivity matrix DOFs

%% DOFs partitioning
DofNrs = reshape(1:n_DOFs, nDOF*(n_ely+1), n_elz+1, n_elx+1);
dofs.L  = DofNrs(nDOF+1:end-nDOF,:,1); dofs.L = dofs.L(:);
dofs.R  = DofNrs(nDOF+1:end-nDOF,:,end); dofs.R = dofs.R(:);
dofs.T  = DofNrs(1:nDOF,:,2:end-1); dofs.T = dofs.T(:);
dofs.B  = DofNrs(end-(nDOF-1):end,:,2:end-1); dofs.B = dofs.B(:);
dofs.TL = DofNrs(1:nDOF,:,1); dofs.TL = dofs.TL(:);
dofs.TR = DofNrs(1:nDOF,:,end); dofs.TR = dofs.TR(:);
dofs.BL = DofNrs(end-(nDOF-1):end,:,1); dofs.BL = dofs.BL(:);
dofs.BR = DofNrs(end-(nDOF-1):end,:,end); dofs.BR = dofs.BR(:);
dofs.I  = DofNrs(nDOF+1:end-nDOF,:,2:end-1); dofs.I = dofs.I(:);

%% System matrix assembly
% Assemble bare plate system matrices
iL = reshape(kron(cMat,ones(24,1))',24*24*n_elem,1);
jL = reshape(kron(cMat,ones(1,24))',24*24*n_elem,1);
sK = reshape(KE(:)*(ones(n_elem,1)'),24*24*n_elem,1);
sM = reshape(ME(:)*(ones(n_elem,1)'),24*24*n_elem,1);
K = sparse(iL,jL,sK); 
M = sparse(iL,jL,sM); 
% Add scatterer (point mass/resonator)
switch scatterer
    case 'mass' 
      mass = m_ratio*mass_UC; 
      dofM = 3*nodeNrs(floor(end/2),end,floor(end/2));
      M(dofM,dofM) = M(dofM,dofM) + mass;
    case 'resonator'
      % Add extra zero row and column to K, M 
      K = [K zeros(n_DOFs,1);zeros(1,n_DOFs) 0];
      M = [M zeros(n_DOFs,1);zeros(1,n_DOFs) 0];
      % Define mass and stiffness
      mass = m_ratio*mass_UC;
      k_res = (f_res*2*pi)^2*mass; 
      dofK = 3*nodeNrs(floor(end/2),end,floor(end/2));
      n_DOFs = n_DOFs+1;
      dofM = n_DOFs;
      % Add resonator matrices to bare plate UC system matrices
      K([dofK dofM],[dofK dofM]) = K([dofK dofM],[dofK dofM]) + [k_res -k_res; -k_res  k_res];
      M([dofK dofM],[dofK dofM]) = M([dofK dofM],[dofK dofM]) + [0 0; 0 mass];
      dofs.I = [dofs.I; dofM];
end

%% Sampling the IBC
mu = 1i*cont_co(1,:).';
tot_steps = zeros(1,size(cont_co,1));
for i = 1 : size(cont_co,1)-1
    n_steps = ceil((sqrt((cont_co(i,1)-cont_co(i+1,1))^2+(cont_co(i,2)-cont_co(i+1,2))^2))/(step_size));
    ed = zeros(2,n_steps);
    for j = 1:2 % j = 1 -> x; j = 2 -> y
        step = (cont_co(i+1,j)-cont_co(i,j))/(n_steps); 
        if step == 0 
            ed(j,:) = cont_co(i,j)*ones(1,n_steps); 
        else 
            ed(j,:) = (cont_co(i,j)+step):step:cont_co(i+1,j); 
        end
    end
    mu = [mu, 1i*ed];
    tot_steps(i+1) = tot_steps(i)+n_steps;
end

%% Dispersion curve calculation
omega = zeros(n_curves,size(mu,2));
for i = 1:size(mu,2)
    % Construction of periodicity matrix
    R = eye(n_DOFs,n_DOFs); 
    R(dofs.R,dofs.L) = exp(mu(1,i))*eye(length(dofs.R),length(dofs.L));
    R(dofs.T,dofs.B) = exp(mu(2,i))*eye(length(dofs.T),length(dofs.B));
    R(dofs.TL,dofs.BL) = exp(mu(2,i))*eye(length(dofs.TL),length(dofs.BL));
    R(dofs.TR,dofs.BL) = exp(mu(1,i)+mu(2,i))*eye(length(dofs.TR),length(dofs.BL));
    R(dofs.BR,dofs.BL) = exp(mu(1,i))*eye(length(dofs.BR),length(dofs.BL));
    R = sparse(R(:,setdiff(1:n_DOFs,[dofs.R;dofs.T;dofs.TL;dofs.BR;dofs.TR])));
    % Impose periodicity boundary conditions
    K_BF = R'*K*R; 
    M_BF = R'*M*R; 
    % Compute dispersion curves
    [~,s] = eigs(K_BF,M_BF,n_curves,0);
    [s,~] = sort(diag(s));
    omega(:,i) = sqrt(s);
end

%% Plot
figure
plot(0:tot_steps(end),real(omega(:,:))/(2*pi),'ko','LineWidth',1);
xlabel('Re(\mu) [-]'); ylabel('Frequency [Hz]'); axis tight;
set(gca,'XTick',tot_steps(:),'XTickLabel',cont_name,'XGrid','on')

%% Function of element matrices
function [KE,ME] = KM(E,nu,rho,a,b,c)
syms x y z
C = (1/((1+nu)*(1-2*nu)))*[1-nu,nu,nu,0,0,0;
    nu,1-nu,nu,0,0,0;
    nu,nu,(1-nu),0,0,0;
    0,0,0,(1-2*nu)/2,0,0;
    0,0,0,0,(1-2*nu)/2,0;
    0,0,0,0,0,(1-2*nu)/2];
N1 = (1/(8*a*b*c))*(a-x)*(b-y)*(c-z);
N2 = (1/(8*a*b*c))*(a+x)*(b-y)*(c-z);
N3 = (1/(8*a*b*c))*(a+x)*(b+y)*(c-z);
N4 = (1/(8*a*b*c))*(a-x)*(b+y)*(c-z);
N5 = (1/(8*a*b*c))*(a-x)*(b-y)*(c+z);
N6 = (1/(8*a*b*c))*(a+x)*(b-y)*(c+z);
N7 = (1/(8*a*b*c))*(a+x)*(b+y)*(c+z);
N8 = (1/(8*a*b*c))*(a-x)*(b+y)*(c+z);

B = [diff(N1,x),0,0,diff(N2,x),0,0,diff(N3,x),0,0,diff(N4,x),0,0,diff(N5,x),0,0,diff(N6,x),0,0,diff(N7,x),0,0,diff(N8,x),0,0;
    0,diff(N1,y),0,0,diff(N2,y),0,0,diff(N3,y),0,0,diff(N4,y),0,0,diff(N5,y),0,0,diff(N6,y),0,0,diff(N7,y),0,0,diff(N8,y),0;
    0,0,diff(N1,z),0,0,diff(N2,z),0,0,diff(N3,z),0,0,diff(N4,z),0,0,diff(N5,z),0,0,diff(N6,z),0,0,diff(N7,z),0,0,diff(N8,z);
    diff(N1,y),diff(N1,x),0,diff(N2,y),diff(N2,x),0,diff(N3,y),diff(N3,x),0,diff(N4,y),diff(N4,x),0,diff(N5,y),diff(N5,x),0,diff(N6,y),diff(N6,x),0,diff(N7,y),diff(N7,x),0,diff(N8,y),diff(N8,x),0;
    0,diff(N1,z),diff(N1,y),0,diff(N2,z),diff(N2,y),0,diff(N3,z),diff(N3,y),0,diff(N4,z),diff(N4,y),0,diff(N5,z),diff(N5,y),0,diff(N6,z),diff(N6,y),0,diff(N7,z),diff(N7,y),0,diff(N8,z),diff(N8,y);
    diff(N1,z),0,diff(N1,x),diff(N2,z),0,diff(N2,x),diff(N3,z),0,diff(N3,x),diff(N4,z),0,diff(N4,x),diff(N5,z),0,diff(N5,x),diff(N6,z),0,diff(N6,x),diff(N7,z),0,diff(N7,x),diff(N8,z),0,diff(N8,x)];
KE = int(int(int(B.'*C*B, z,-c,c),y,-b,b),x,-a,a);
KE = E*double(KE);
[~, kua, kaa, kau] = FEM_incompatible_modes(E,nu,a,b,c);
KE = KE - kua*(kaa\kau);
KE = KE/2;

N = [N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8,0,0;
     0,N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8,0;
     0,0,N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8];
ME = int(int(int(N.'*N, z,-c,c),y,-b,b),x,-a,a);
ME = rho*double(ME/8);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This Matlab code was written in October 2023                        %%%
%%% by V. Cool, E. Deckers, L. Van Belle and C. Claeys                  %%%
%%% Department of Mechanical Engineering, 3001 Heverlee, Belgium        %%%
%%% Any comments are welcome: claus.claeys@kuleuven.be                  %%%
%%%                                                                     %%%
%%% The code is intended for educational purposes.                      %%%
%%% Theoretical details are discussed in the corresponding paper:       %%%
%%% "A guide to numerical dispersion curve calculations:                %%%
%%% explanation, interpretation and basic Matlab code"                  %%%
%%%                                                                     %%%
%%% This code can be dowloaded from the corresponding Github page:      %%%
%%% github.com/LMSD-KULeuven/2D_InverseUndamped_DispersionCurves        %%%
%%%                                                                     %%%
%%%                                                                     %%%
%%% Copyright (c) 2023 KU Leuven Mecha(tro)nic System Dynamics (LMSD)   %%%
%%%                                                                     %%%
%%% Permission is hereby granted, free of charge, to any person         %%%
%%% obtaining a copy of this software and associated documentation      %%%
%%% files (the "Software"), to deal in the Software without             %%%
%%% restriction, including without limitation the rights to use,        %%%
%%% copy, modify, merge, publish, distribute, sublicense, and/or sell   %%%
%%% copies of the Software, and to permit persons to whom the           %%%
%%% Software is furnished to do so, subject to the following            %%%
%%% conditions:                                                         %%%
%%%                                                                     %%%
%%% The above copyright notice and this permission notice shall be      %%%
%%% included in all copies or substantial portions of the Software.     %%%
%%%                                                                     %%%
%%% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,     %%%
%%% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES     %%%
%%% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND            %%%
%%% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT         %%%
%%% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,        %%%
%%% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING        %%%
%%% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR       %%%
%%% OTHER DEALINGS IN THE SOFTWARE.                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


