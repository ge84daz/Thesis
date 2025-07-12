%__________________________________________________________________________
%
%        f_planequad.m    
%
%       INPUT: 
%       ey = [y1 y2 y3]         element coordinates
%       ez = [z1 z2 z3]
%       ep = [ptype t ]         ptype: analysis type
%                               t: thickness
%       D                       constitutive matrix
%
%       eq = [bx;               bx: body force x-dir
%             by]               by: body force y-dir
%       
%       OUTPUT: 
%       Ke : element stiffness matrix (6 x 6)
%       fe : equivalent nodal forces (6 x 1)
%
%       DESCRIPTION: 
%       The function calculates the stiffness matrix of a four-node plane
%       stress or strain element within the image space
%
%       REMARK: 
%       Original Author: Manuela Hackenberg (planqe.m)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

  
function Ke = f_planequad(ey,ez,ep,ngp,wgaus,zeta,eta,Kx,Omega,D1,rho1)

ptype   = ep(1);
d       = ep(2);

bx      = 0; 
by      = 0.; 
if nargin==10  
    bx  =   eq(1); 
    by   =  eq(2); 
end

N=zeros(ngp); % shape functions - Diss Hackenberg eq. 4.11
N(:,1)=(1-zeta).*(1-eta)/4;  N(:,2)=(1+zeta).*(1-eta)/4;
N(:,3)=(1+zeta).*(1+eta)/4;  N(:,4)=(1-zeta).*(1+eta)/4;
NMe=zeros(12,12);

r2=ngp*2;
dNr=zeros(2*ngp,4);
dNr(1:2:r2,1)=-(1-eta)/4;       dNr(1:2:r2,2)=+(1-eta)/4;       dNr(1:2:r2,3)= +(1+eta)/4;      dNr(1:2:r2,4)= -(1+eta)/4;
dNr(2:2:r2,1)=-(1-zeta)/4;      dNr(2:2:r2,2)= -(1+zeta)/4;     dNr(2:2:r2,3)= +(1+zeta)/4;     dNr(2:2:r2,4)=+(1-zeta)/4;

J=dNr*[ey;ez]';
Ke=complex(zeros(12,12));
Me=complex(zeros(12,12));

for s=1:ngp
    detJ=det(J((s*2-1):s*2,:));
    dNyz=J((s*2-1):s*2,:)\dNr((s*2-1):s*2,:);
    dNx=complex(zeros(1,ngp)); % predefinition
    dNx=1i*Kx*N(s,:);
    
    B=complex(zeros(6,12));
    for t =1:4
    B(1,t*3-2)=dNx(t);
    B(2,t*3-1)=dNyz(1,t);
    B(3,t*3)=dNyz(2,t);
    B(4,t*3-2)=dNyz(1,t);
    B(4,t*3-1)=dNx(t);
    B(5,t*3-1)=dNyz(2,t);
    B(5,t*3)=dNyz(1,t);
    B(6,t*3-2)=dNyz(2,t);
    B(6,t*3)=dNx(t);
    end
    
    Ke=Ke+B'*D1*B*detJ*wgaus(s)*d;
    NMe((s*3-2),1:3:3*ngp)=N(s,:);
    NMe((s*3-1),2:3:3*ngp)=N(s,:);
    NMe((s*3-0),3:3:3*ngp)=N(s,:);
    Me=Me+NMe((s*3-2):s*3,:)'*rho1*NMe((s*3-2):s*3,:)*detJ*wgaus(s)*d;
    
  
end 

%% Output
Ke=Ke-Omega^2*Me;


end
