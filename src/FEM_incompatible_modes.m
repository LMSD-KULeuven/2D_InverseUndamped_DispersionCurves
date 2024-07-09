%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%             INCOMPATIBLE MODE ELEMENT SYSTEM MATRICES            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kuu, kua, kaa, kau] = FEM_incompatible_modes(E,nux,ax,bx,cx)
% INPUTS
% E                      Young's modulus [Pa]
% nux                    Poisson's ratio [-]
% ax, bx, cx             The size of the element in the x,y and z-direction.

ncoord = 3;              % number of coordinates: 3 for 3D
ndof = 3;                % number of degrees-of-freedom per node (often ndof = ncoord)
nelnodes = 8;            % number of nodes per element
npoints = 8;             % number of integration points: hard coded for 3D solids
G = E/(2*(1+nux));
materialprops = [G;nux]; % material properties: constitutive procedures
coord = [-ax -bx -cx;
          ax -bx -cx;
          ax  bx -cx;
         -ax  bx -cx;
         -ax -bx  cx;
          ax -bx  cx;
          ax  bx  cx; 
         -ax  bx  cx]';  % coords(i,a): ith coord of ath node  

% Initialization element matrices
kuu = zeros(ndof*nelnodes,ndof*nelnodes);
kau = zeros(ncoord*ndof,ndof*nelnodes);
kua = zeros(ndof*nelnodes,ncoord*ndof);
kaa = zeros(ncoord*ndof,ncoord*ndof);
     
% Local variables          
% xi(i,inpt)         ith local coord of integration point no. intpt
% N(a)               Shape function associated with ath node on element
% dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
% dNdx(a,i)          Derivative of ath shape function wrt ith global coord
% dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
% dxidx(i,j)         Derivative of ith local coord wrt jth global coord
% det                Determinant of jacobian

% Integration weights
w = ones(1, npoints);

% Integration points
x1D = [-0.5773502692,0.5773502692];
[x1, x2, x3] = ndgrid(x1D, x1D, x1D);
xilist = [x1(:)'; x2(:)'; x3(:)'];

% Compute jacobian at center of element
xi = zeros(ncoord, 1);
dNdxi = shapefunctionderivs(xi);
dxdxi = coord * dNdxi;
dt0 = det(dxdxi);

% Compute the elasticity matrix
mu = materialprops(1);
nu = materialprops(2);
dsde = ((mu*(2*(1+nu)))/((1+nu)*(1-2*nu)))*...
    [1-nu,nu,nu,0,0,0;
    nu,1-nu,nu,0,0,0;
    nu,nu,1-nu,0,0,0;
    0,0,0,(1-2*nu)/2,0,0;
    0,0,0,0,(1-2*nu)/2,0;
    0,0,0,0,0,(1-2*nu)/2];

% Loop over the integration points
for intpt = 1:npoints
    % Local coordinate integration point
    xi = xilist(:, intpt);
    
    % Compute the jacobian matrix & its determinant
    dNdxi = shapefunctionderivs(xi);
    dxdxi = coord * dNdxi;
    dxidx = inv(dxdxi);
    dt = det(dxdxi);
    
    % Convert shape function derivatives: derivatives wrt global coords
    dNdx = dNdxi * dxidx;
    
    % Combination
    dNdx_full = [dNdx(1,1) 0 0 dNdx(2,1) 0 0 dNdx(3,1) 0 0 dNdx(4,1) 0 0 dNdx(5,1) 0 0 dNdx(6,1) 0 0 dNdx(7,1) 0 0 dNdx(8,1) 0 0;
        0 dNdx(1,2) 0 0 dNdx(2,2) 0 0 dNdx(3,2) 0 0 dNdx(4,2) 0 0 dNdx(5,2) 0 0 dNdx(6,2) 0 0 dNdx(7,2) 0 0 dNdx(8,2) 0;
        0 0 dNdx(1,3) 0 0 dNdx(2,3) 0 0 dNdx(3,3) 0 0 dNdx(4,3) 0 0 dNdx(5,3) 0 0 dNdx(6,3) 0 0 dNdx(7,3) 0 0 dNdx(8,3);
        dNdx(1,2) dNdx(1,1) 0 dNdx(2,2) dNdx(2,1) 0 dNdx(3,2) dNdx(3,1) 0 dNdx(4,2) dNdx(4,1) 0 dNdx(5,2) dNdx(5,1) 0 dNdx(6,2) dNdx(6,1) 0 dNdx(7,2) dNdx(7,1) 0 dNdx(8,2) dNdx(8,1) 0;
        0  dNdx(1,3) dNdx(1,2) 0 dNdx(2,3) dNdx(2,2) 0 dNdx(3,3) dNdx(3,2) 0 dNdx(4,3) dNdx(4,2) 0 dNdx(5,3) dNdx(5,2) 0 dNdx(6,3) dNdx(6,2) 0 dNdx(7,3) dNdx(7,2) 0 dNdx(8,3) dNdx(8,2);
        dNdx(1,3) 0 dNdx(1,1) dNdx(2,3) 0 dNdx(2,1) dNdx(3,3) 0 dNdx(3,1) dNdx(4,3) 0 dNdx(4,1) dNdx(5,3) 0 dNdx(5,1) dNdx(6,3) 0 dNdx(6,1) dNdx(7,3) 0 dNdx(7,1) dNdx(8,3) 0 dNdx(8,1)];
    
    % Correction matrix
    dBdx_full = (dt0/dt)*[xi(1)*dxidx(1,1) 0 0 xi(2)*dxidx(2,1) 0 0 xi(3)*dxidx(3,1) 0 0;
        0 xi(1)*dxidx(1,2) 0 0 xi(2)*dxidx(2,2) 0 0 xi(3)*dxidx(3,2) 0;
        0 0 xi(1)*dxidx(1,3) 0 0 xi(2)*dxidx(2,3) 0 0 xi(3)*dxidx(3,3);
        xi(1)*dxidx(1,2) xi(1)*dxidx(1,1) 0 xi(2)*dxidx(2,2) xi(2)*dxidx(2,1) 0 xi(3)*dxidx(3,2) xi(3)*dxidx(3,1) 0;
        0  xi(1)*dxidx(1,3) xi(1)*dxidx(1,2) 0 xi(2)*dxidx(2,3) xi(2)*dxidx(2,2) 0 xi(3)*dxidx(3,3) xi(3)*dxidx(3,2);
        xi(1)*dxidx(1,3) 0 xi(1)*dxidx(1,1) xi(2)*dxidx(2,3) 0 xi(2)*dxidx(2,1) xi(3)*dxidx(3,3) 0 xi(3)*dxidx(3,1)];
    
    % Element matrices
    kuu = kuu + dNdx_full.'*dsde*dNdx_full*w(intpt)*dt;
    kau = kau + dBdx_full.'*dsde*dNdx_full*w(intpt)*dt;
    kua = kua + dNdx_full.'*dsde*dBdx_full*w(intpt)*dt;
    kaa = kaa + dBdx_full.'*dsde*dBdx_full*w(intpt)*dt;

end  
end

function dNdxi = shapefunctionderivs(xi)
dNdxi = zeros(8,3);
dNdxi(1,1) = -(1.-xi(2))*(1.-xi(3))/8.; %dN1 dx1
dNdxi(1,2) = -(1.-xi(1))*(1.-xi(3))/8.; %dN1 dx2
dNdxi(1,3) = -(1.-xi(1))*(1.-xi(2))/8.; %dN1 dx3
dNdxi(2,1) =  (1.-xi(2))*(1.-xi(3))/8.; %dN2 dx1
dNdxi(2,2) = -(1.+xi(1))*(1.-xi(3))/8.; %dN2 dx2
dNdxi(2,3) = -(1.+xi(1))*(1.-xi(2))/8.; %dN2 dx3
dNdxi(3,1) =  (1.+xi(2))*(1.-xi(3))/8.; %...
dNdxi(3,2) =  (1.+xi(1))*(1.-xi(3))/8.;
dNdxi(3,3) = -(1.+xi(1))*(1.+xi(2))/8.;
dNdxi(4,1) = -(1.+xi(2))*(1.-xi(3))/8.;
dNdxi(4,2) =  (1.-xi(1))*(1.-xi(3))/8.;
dNdxi(4,3) = -(1.-xi(1))*(1.+xi(2))/8.;
dNdxi(5,1) = -(1.-xi(2))*(1.+xi(3))/8.;
dNdxi(5,2) = -(1.-xi(1))*(1.+xi(3))/8.;
dNdxi(5,3) =  (1.-xi(1))*(1.-xi(2))/8.;
dNdxi(6,1) =  (1.-xi(2))*(1.+xi(3))/8.;
dNdxi(6,2) = -(1.+xi(1))*(1.+xi(3))/8.;
dNdxi(6,3) =  (1.+xi(1))*(1.-xi(2))/8.;
dNdxi(7,1) =  (1.+xi(2))*(1.+xi(3))/8.;
dNdxi(7,2) =  (1.+xi(1))*(1.+xi(3))/8.;
dNdxi(7,3) =  (1.+xi(1))*(1.+xi(2))/8.;
dNdxi(8,1) = -(1.+xi(2))*(1.+xi(3))/8.;
dNdxi(8,2) =  (1.-xi(1))*(1.+xi(3))/8.;
dNdxi(8,3) =  (1.-xi(1))*(1.+xi(2))/8.;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This Matlab code was written in July 2024                        
%%% by V. Cool, E. Deckers, L. Van Belle and C. Claeys                  
%%% Department of Mechanical Engineering, 3001 Heverlee, Belgium        
%%% Any comments are welcome: claus.claeys@kuleuven.be                  
%%%                                                                     
%%% The code is inspired on the implementation found on       
%%% https://solidmechanics.org/FEA.php                                              
%%%                                                                     
%%% This code can be dowloaded from the corresponding Github page:      
%%% github.com/LMSD-KULeuven/2D_InverseUndamped_DispersionCurves             
%%%                                                                                                                                         
%%% Copyright (c) 2024 KU Leuven Mecha(tro)nic System Dynamics (LMSD)        
%%%                                                                     
%%% Permission is hereby granted, free of charge, to any person         
%%% obtaining a copy of this software and associated documentation      
%%% files (the "Software"), to deal in the Software without             
%%% restriction, including without limitation the rights to use,        
%%% copy, modify, merge, publish, distribute, sublicense, and/or sell   
%%% copies of the Software, and to permit persons to whom the           
%%% Software is furnished to do so, subject to the following            
%%% conditions:                                                         
%%%                                                                     
%%% The above copyright notice and this permission notice shall be      
%%% included in all copies or substantial portions of the Software.     
%%%                                                                     
%%% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,     
%%% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES     
%%% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND            
%%% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT         
%%% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,        
%%% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING        
%%% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR       
%%% OTHER DEALINGS IN THE SOFTWARE.                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

