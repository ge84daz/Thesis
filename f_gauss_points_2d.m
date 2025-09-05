%__________________________________________________________________________
%
%        f_gauss_points_2D.m    
%
%       INPUT: 
%       - nPoints               = Number of Gauss Points for each direction
%
%       OUTPUT: 
%       - wp(nPoints^2)         = Weights of Points
%       - gp(nPoints^2,2)       = Coordinates of Points
%
%       DESCRIPTION: 
%       - Set up Gauss points (weights and coordinates)
%
%       REMARK: 
%       Original Author: Manuela Hackenberg (gauss_points_2D)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [wp,xsi,eta,ngp] = f_gauss_points_2d (nPoints)

if nPoints==1
    g1=0.0; w1=2.0;
    xsi= g1;
    eta = g1;
    w=zeros(nPoints*1,2);
    w=[ w1 w1 ];
    
elseif nPoints==2 % Diss Hackenberg p.58
    g1=0.577350269189626; w1=1;
    xsi=[-g1; g1;-g1; g1];
    eta=[-g1;-g1; g1; g1];
    w=zeros(nPoints*2,2);
    w(:,1)=[ w1; w1; w1; w1];
    w(:,2)=[ w1; w1; w1; w1];
    
elseif nPoints==3
    g1=0.774596669241483;
    g2=0.0;
    w1=0.555555555555555;
    w2=0.888888888888888;
    xsi=[-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    eta=[-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    w=zeros(3*nPoints,2);
    w(:,1)=[ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2)=[ w1; w1; w1; w2; w2; w2; w1; w1; w1];
    
else
    disp('Used number of integration points not implemented');
    xsi=0;
    eta=0;
    ngp=0;
    wp=[0 0];
    return
    
end


wp=w(:,1).*w(:,2);
ngp=nPoints^2;
end % end function





