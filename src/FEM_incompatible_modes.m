%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This Matlab code was dowloaded from                                 %%%
%%% https://solidmechanics.org/FEA.php                                  %%% 
%%% in October 2023 and adapted to compute the incompatible mode        %%%
%%% element system matrices of 3D linear solid elements.                %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kuu, kua, kaa, kau] = FEM_incompatible_modes(E,nux,ax,bx,cx)
%%%%% INPUTS
% E                   Young's modulus [Pa]
% nux                 Poisson's ratio [-]
% ax, bx, cx          The size of the element in the x,y and z-direction.


ncoord = 3;             % 3 for 3D
ndof = 3;               % number of degrees-of-freedom per node
nelnodes = 8;           % number of nodes per element
G = E/(2*(1+nux));
materialprops = [G;nux];
coord = [-ax,-bx,-cx;
          ax -bx -cx;
          ax  bx -cx;
         -ax  bx -cx;
         -ax -bx  cx;
          ax -bx  cx;
          ax  bx  cx; 
         -ax  bx  cx]';
[kuu, kua, kaa, kau] = elstif(ncoord,ndof,nelnodes,coord,materialprops);

function [kuu, kua, kaa, kau] = elstif(ncoord,ndof,nelnodes,coord,materialprops)
%================= ELEMENT STIFFNESS MATRIX ================================
%  Assemble the element stiffness
%   Arguments;
%      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
%      ndof               No. degrees of freedom per node (often ndof = ncoord)
%      nelnodes           No. nodes on the element
%      coords(i,a)        ith coord of ath node
%      materialprops      Material properties passed on:constitutive procedures
%      displacement(i,a)  ith displacement component at ath node
%
%   Local variables
%      npoints            No. integration points
%      xi(i,inpt)         ith local coord of integration point no. intpt
%      w(intpt)           weight for integration point no. intpt
%      N(a)               Shape function associated with ath node on element
%      dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
%      dNdx(a,i)          Derivative of ath shape function wrt ith global coord
%      dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
%      dxidx(i,j)         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian
%      strain(i,j)        strain_ij components
%      dsde(i,j,k,l)      Derivative of stress_ij with respect:strain_kl
%      kel(row,col)       Rows && cols of element stiffness

   npoints = 8; % Hard Coded for 3D solids
   dNdx = zeros(nelnodes,ncoord);
   dxdxi = zeros(ncoord,ncoord);
   kuu = zeros(ndof*nelnodes,ndof*nelnodes);
   kau = zeros(ncoord*ndof,ndof*nelnodes);
   kua = zeros(ndof*nelnodes,ncoord*ndof);
   kaa = zeros(ncoord*ndof,ncoord*ndof);

   % Integration weights
    w = [1.,1.,1.,1.,1.,1.,1.,1.];
   % Integration points
   xilist = zeros(ncoord,npoints);
   x1D = [-0.5773502692,0.5773502692];
   for k = 1:2
       for j = 1:2
           for i = 1:2
               n = 4*(k-1) + 2*(j-1) + i;
               xilist(1,n) = x1D(i);
               xilist(2,n) = x1D(j);
               xilist(3,n) = x1D(k);
           end
       end
   end
  
   % Compute jacobian at center of element
   xi(1:ncoord) = 0;
   dNdxi = shapefunctionderivs(xi);
   for i = 1:ncoord
     for j = 1:ncoord
       dxdxi(i,j) = 0.;
       for a = 1:nelnodes
         dxdxi(i,j) = dxdxi(i,j) + coord(i,a)*dNdxi(a,j);
       end
     end
   end
   dt0 = det(dxdxi);
   
   % Loop over the integration points
   for intpt = 1:npoints
      % Compute shape functions && derivatives wrt local coords
      for i = 1:ncoord
        xi(i) = xilist(i,intpt);
      end      
      dNdxi = shapefunctionderivs(xi);
      % Compute the jacobian matrix && its determinant
      for i = 1:ncoord
        for j = 1:ncoord
          dxdxi(i,j) = 0.;
          for a = 1:nelnodes
            dxdxi(i,j) = dxdxi(i,j) + coord(i,a)*dNdxi(a,j);
          end
        end
      end
      dxidx = inv(dxdxi);
      dt = det(dxdxi);
      % Convert shape function derivatives:derivatives wrt global coords
      for a = 1:nelnodes
        for i = 1:ncoord
          dNdx(a,i) = 0.;
          for j = 1:ncoord
            dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
          end
        end
      end
      % Compute the material tangent stiffness (d stress/d strain)
      % ds/de is just C_ijkl for linear elasticity
      dsde = materialstiffness(ndof,ncoord,materialprops);
      % Compute the element stiffness        
      for a = 1:nelnodes
        for i = 1:ndof
          for b = 1:nelnodes
            for k = 1:ndof
              row = ndof*(a-1)+i;
              col = ndof*(b-1)+k;
              for j = 1:ncoord
                for l = 1:ncoord
                  kuu(col,row) = kuu(col,row) + dsde(i,j,k,l)*dNdx(b,l)*dNdx(a,j)*w(intpt)*dt;
                end
              end
            end
          end
          
          for m = 1 : ncoord
            for k = 1 : ncoord
              row = ndof*(a-1)+i;
              col = (m-1)*ncoord + k;
              for j = 1 : ncoord
                for l = 1 : ncoord
                   kau(col,row) = kau(col,row) + ...
                                 dsde(i,j,k,l)*(xi(m)*dt0/dt)*dxidx(m,l)*dNdx(a,j)*w(intpt)*dt;
                   kua(row,col) = kua(row,col) + ...
                                 dsde(i,j,k,l)*(xi(m)*dt0/dt)*dxidx(m,l)*dNdx(a,j)*w(intpt)*dt;
                end
              end
            end
          end
        end
      end
      for n = 1 : ndof
        for i = 1 : ncoord
          for m = 1 : ndof
            for k = 1 : ncoord
              row = (n-1)*ncoord + i;
              col = (m-1)*ncoord + k;
              for j = 1 : ncoord
                for l = 1 : ncoord
                   kaa(col,row) = kaa(col,row) + ...
                              dsde(i,j,k,l)*(xi(m)*dt0/dt)*dxidx(m,l)*(xi(n)*dt0/dt)*dxidx(n,j)*w(intpt)*dt;
                end
              end
            end
          end        
        end
      end
   end
end

function C = materialstiffness(ndof,ncoord,materialprops)
   % Computes elasticity tensor C_{ijkl} = shear modulus and Poissons ratio
   mu = materialprops(1);
   nu = materialprops(2);
   C = zeros(ndof,ncoord,ndof,ncoord);
   if (ncoord == 3) 
     for i = 1:3 
       for j = 1:3
         for k = 1:3
           for l = 1:3
             if (i==j && k==l)
                 C(i,j,k,l)=C(i,j,k,l) + 2.*mu*nu/(1.-2.*nu);
             end
             if (i==k && j==l)
                 C(i,j,k,l)=C(i,j,k,l) + mu;
             end
             if (i==l && j==k)
                 C(i,j,k,l)=C(i,j,k,l) + mu;
             end
            end
          end
        end
      end
   end
end

function dNdxi = shapefunctionderivs(xi)
    dNdxi = zeros(8,3);
    dNdxi(1,1) = -(1.-xi(2))*(1.-xi(3))/8.;
    dNdxi(1,2) = -(1.-xi(1))*(1.-xi(3))/8.;
    dNdxi(1,3) = -(1.-xi(1))*(1.-xi(2))/8.;
    dNdxi(2,1) = (1.-xi(2))*(1.-xi(3))/8.;
    dNdxi(2,2) = -(1.+xi(1))*(1.-xi(3))/8.;
    dNdxi(2,3) = -(1.+xi(1))*(1.-xi(2))/8.;
    dNdxi(3,1) = (1.+xi(2))*(1.-xi(3))/8.;
    dNdxi(3,2) = (1.+xi(1))*(1.-xi(3))/8.;
    dNdxi(3,3) = -(1.+xi(1))*(1.+xi(2))/8.;
    dNdxi(4,1) = -(1.+xi(2))*(1.-xi(3))/8.;
    dNdxi(4,2) = (1.-xi(1))*(1.-xi(3))/8.;
    dNdxi(4,3) = -(1.-xi(1))*(1.+xi(2))/8.;
    dNdxi(5,1) = -(1.-xi(2))*(1.+xi(3))/8.;
    dNdxi(5,2) = -(1.-xi(1))*(1.+xi(3))/8.;
    dNdxi(5,3) = (1.-xi(1))*(1.-xi(2))/8.;
    dNdxi(6,1) = (1.-xi(2))*(1.+xi(3))/8.;
    dNdxi(6,2) = -(1.+xi(1))*(1.+xi(3))/8.;
    dNdxi(6,3) = (1.+xi(1))*(1.-xi(2))/8.;
    dNdxi(7,1) = (1.+xi(2))*(1.+xi(3))/8.;
    dNdxi(7,2) = (1.+xi(1))*(1.+xi(3))/8.;
    dNdxi(7,3) = (1.+xi(1))*(1.+xi(2))/8.;
    dNdxi(8,1) = -(1.+xi(2))*(1.+xi(3))/8.;
    dNdxi(8,2) = (1.-xi(1))*(1.+xi(3))/8.;
    dNdxi(8,3) = (1.-xi(1))*(1.+xi(2))/8.;
end
end