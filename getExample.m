
function ex = getExample(example,n)
ex = zeros(n,n);
%-------------------------------------------------------------------------
% One source and two sinks 
if ( strcmp(example,'3Points') )
    M = [round(n/2),round(n/2)];
    r = round(n/2)-round(n/6);
    P1 = round(r*[sin(0),cos(0)])+M;
    P2 = round(r*[sin(2*pi/3),cos(2*pi/3)])+M;
    P3 = round(r*[sin(4*pi/3),cos(4*pi/3)])+M;
    ex(P1(1),P1(2)) = 1*(n-1)^2;  
    ex(P2(1),P2(2)) = -1/2*(n-1)^2;
    ex(P3(1),P3(2)) = -1/2*(n-1)^2;
%-------------------------------------------------------------------------
% One source and three sinks 
elseif ( strcmp(example,'4Points') ) 
    M = [round(n/2),round(n/2)];
    r = round(n/2)-round(n/6);
    P1 = round(r*[sin(pi/4),cos(pi/4)])+M+[-3,3];
    P2 = round(r*[sin(3*pi/4),cos(3*pi/4)])+M;
    P3 = round(r*[sin(5*pi/4),cos(5*pi/4)])+M;
    P4 = round(r*[sin(7*pi/4),cos(7*pi/4)])+M;
    ex(P1(1),P1(2)) = 1*(n-1)^2;
    ex(P2(1),P2(2)) = -1/3*(n-1)^2;
    ex(P3(1),P3(2)) = -1/3*(n-1)^2;   
    ex(P4(1),P4(2)) = -1/3*(n-1)^2; 
%-------------------------------------------------------------------------  
% One source and four sinks 
elseif ( strcmp(example,'5Points') )
    M = [round(n/2),round(n/2)];
    r = round(n/2)-round(n/6);
    P1 = round(r*[sin(0),cos(0)])+M;
    P2 = round(r*[sin(2*pi/5),cos(2*pi/5)])+M+[-3,3];
    P3 = round(r*[sin(4*pi/5),cos(4*pi/5)])+M;
    P4 = round(r*[sin(6*pi/5),cos(6*pi/5)])+M;
    P5 = round(r*[sin(8*pi/5),cos(8*pi/5)])+M;
    ex(P1(1),P1(2)) = 1*(n-1)^2;
    ex(P2(1),P2(2)) = -1/4*(n-1)^2;
    ex(P3(1),P3(2)) = -1/4*(n-1)^2;   
    ex(P4(1),P4(2)) = -1/4*(n-1)^2;
    ex(P5(1),P5(2)) = -1/4*(n-1)^2;  
%-------------------------------------------------------------------------  
% One source and five sinks   
elseif ( strcmp(example,'6Points') )
    M = [round(n/2),round(n/2)];
    r = round(n/2)-round(n/6);
    P1 = round(r*[sin(0),cos(0)])+M;
    P2 = round(r*[sin(2*pi/6),cos(2*pi/6)])+M;
    P3 = round(r*[sin(4*pi/6),cos(4*pi/6)])+M;
    P4 = round(r*[sin(6*pi/6),cos(6*pi/6)])+M;
    P5 = round(r*[sin(8*pi/6),cos(8*pi/6)])+M;
    P6 = round(r*[sin(10*pi/6),cos(10*pi/6)])+M;
    ex(P1(1),P1(2)) = -1/5*(n-1)^2;
    ex(P2(1),P2(2)) = -1/5*(n-1)^2;
    ex(P3(1),P3(2)) = -1/5*(n-1)^2;   
    ex(P4(1),P4(2)) = -1/5*(n-1)^2;
    ex(P5(1),P5(2)) = -1/5*(n-1)^2; 
    ex(P6(1),P6(2)) = 1*(n-1)^2;
%-------------------------------------------------------------------------  
% Four sources and four sinks 
elseif ( strcmp(example,'4To4') )    
    ex(round(n/5),round(n/6)) = 1/4*(n-1)^2;
    ex(round(2*n/5),round(n/6)) = 1/4*(n-1)^2;
    ex(round(3*n/5),round(n/6)) = 1/4*(n-1)^2;
    ex(round(4*n/5),round(n/6)) = 1/4*(n-1)^2;
    ex(round(n/5),round(5*n/6)) = -1/4*(n-1)^2;
    ex(round(2*n/5),round(5*n/6)) = -1/4*(n-1)^2;
    ex(round(3*n/5),round(5*n/6)) = -1/4*(n-1)^2;
    ex(round(4*n/5),round(5*n/6)) = -1/4*(n-1)^2;
%-------------------------------------------------------------------------  
% One source in the middle 
elseif ( strcmp(example,'1To16Circle') )
    numSinks = 16;
    M = [round(n/2),round(n/2)];
    r = round(n/2)-round(n/8);
    ex(round(n/2),round(n/2)) = (n-1)^2;
    P = zeros(numSinks,2);
    x1 = 3/8;
    tau = (2*pi-6)/4;
    P(1,:) = round(r*[sin(x1/2),cos(x1/2)])+M;
    P(2,:) = round(r*[sin(x1/2+x1),cos(x1/2+x1)])+M;
    P(3,:) = round(r*[sin(x1/2+2*x1+tau),cos(x1/2+2*x1+tau)])+M;
    P(4,:) = round(r*[sin(x1/2+3*x1+tau),cos(x1/2+3*x1+tau)])+M;
    P(5,:) = round(r*[sin(x1/2+4*x1+tau),cos(x1/2+4*x1+tau)])+M;
    P(6,:) = round(r*[sin(x1/2+5*x1+tau),cos(x1/2+5*x1+tau)])+M;
    P(7,:) = round(r*[sin(x1/2+6*x1+2*tau),cos(x1/2+6*x1+2*tau)])+M;
    P(8,:) = round(r*[sin(x1/2+7*x1+2*tau),cos(x1/2+7*x1+2*tau)])+M;
    P(9,:) = round(r*[sin(x1/2+8*x1+2*tau),cos(x1/2+8*x1+2*tau)])+M;
    P(10,:) = round(r*[sin(x1/2+9*x1+2*tau),cos(x1/2+9*x1+2*tau)])+M;
    P(11,:) = round(r*[sin(x1/2+10*x1+3*tau),cos(x1/2+10*x1+3*tau)])+M;
    P(12,:) = round(r*[sin(x1/2+11*x1+3*tau),cos(x1/2+11*x1+3*tau)])+M;
    P(13,:) = round(r*[sin(x1/2+12*x1+3*tau),cos(x1/2+12*x1+3*tau)])+M;
    P(14,:) = round(r*[sin(x1/2+13*x1+3*tau),cos(x1/2+13*x1+3*tau)])+M;
    P(15,:) = round(r*[sin(x1/2+14*x1+4*tau),cos(x1/2+14*x1+4*tau)])+M;
    P(16,:) = round(r*[sin(x1/2+15*x1+4*tau),cos(x1/2+15*x1+4*tau)])+M;
    for k = 1:numSinks
        ex(P(k,1),P(k,2)) = -1/numSinks*(n-1)^2;
    end
%-------------------------------------------------------------------------     
% Source in the middle and diffuse sinks outside of circle   
elseif ( strcmp(example,'CircleDiffuse') )
    M = [round(n/2),round(n/2)];
    r = round(n/2)-round(n/10);
    for i = 1:n
        for j = 1:n
            if ( sqrt((i-M(1))^2+(j-M(2))^2) > r )
                ex(i,j) = -1;
            end
        end
    end
    numSinks = sum(abs(ex(:)));
    ex = ex/numSinks*(n-1)^2;
    ex(round(n/2),round(n/2)) = (n-1)^2;
%-------------------------------------------------------------------------  
else
    disp('Example not implemented');
end
%-------------------------------------------------------------------------
end