% Evaluating nonlinear impairments in the optical-fiber channel by implementing the EGN model 

function [gn_chi1,chi11,chi21,chi13,chi23,chi14,chi24,chi34,NLIN_var...
    ,NLIN_var_12,NLIN_var_3,NLIN_var_4] = ...
    calc_interChannel(n_ch,gamma,beta2,alpha,Nspan,L,PD,P0,kur,kur3,N,q)%,
R = 2*pi*(rand(5, N)-0.5*ones(5, N));

Volume = (2*pi)^4;

% Calculate chi1
gn_chi1=zeros(n_ch,n_ch,n_ch,3);
chi11=zeros(n_ch,n_ch,n_ch,3);
chi21=zeros(n_ch,n_ch,n_ch,3);
chi13=zeros(n_ch,n_ch,n_ch,3);
chi23=zeros(n_ch,n_ch,n_ch,3);
chi14=zeros(n_ch,n_ch,n_ch,3);
chi24=zeros(n_ch,n_ch,n_ch,3);
chi34=zeros(n_ch,n_ch,n_ch,3);

% GN_term
for l=-1:1
    for n1=1:n_ch
        for n2=n1:n_ch
            if(1<=n1+n2-1+l && n1+n2-1+l<=n_ch && n1-1+l~=0 && n2-1+l~=0 && n1~=n2)
                w0 = R(1,:)+R(2,:)-R(3,:)+l*2*pi*q;%*l*50e9=*2*pi*q
                arg1 = (R(1,:)-R(3,:)+((n1-1)*2*pi*q)).*(R(2,:)-R(3,:)+((n2-1)*2*pi*q));%*50e9=2*pi*q
                argPD1 = arg1;
                ss1 = exp(1i*argPD1*PD).*(exp(1i*beta2*arg1*L-alpha*L)-1)...
                    ./(1i*beta2*arg1-alpha).*(w0<pi).*(w0>-pi);%-16e9=pi
                s1 = abs(ss1.*(1-exp(1i*Nspan*arg1*beta2*L))...
                    ./(1-exp(1i*arg1*beta2*L))).^2/Volume;
                avgF1 = sum(s1)/N;
                gn_chi1(n1,n2,1,l+2) = (16/27)*avgF1*Volume*(gamma^2*P0^3);
                gn_chi1(n2,n1,1,l+2)=gn_chi1(n1,n2,1,l+2);
                
            end
        end
    end
    
end

for l=-1:1
    for n=2:n_ch
        gn_chi1(2:n_ch,2:n_ch,n,l+2)=gn_chi1(1:n_ch-1,1:n_ch-1,n-1,l+2);
        for n2=1:n_ch
            if(1<=1+n2-n+l && 1+n2-n+l<=n_ch && 1-n+l~=0 && n2-n+l~=0 && 1~=n2)
                w0 = R(1,:)+R(2,:)-R(3,:)+l*2*pi*q;%*l*50e9=*2*pi*q
                arg1 = (R(1,:)-R(3,:)+((1-n)*2*pi*q)).*(R(2,:)-R(3,:)+((n2-n)*2*pi*q));%*50e9=2*pi*q
                argPD1 = arg1;
                ss1 = exp(1i*argPD1*PD).*(exp(1i*beta2*arg1*L-alpha*L)-1)...
                    ./(1i*beta2*arg1-alpha).*(w0<pi).*(w0>-pi);%-16e9=pi
                s1 = abs(ss1.*(1-exp(1i*Nspan*arg1*beta2*L))...
                    ./(1-exp(1i*arg1*beta2*L))).^2/Volume;
                avgF1 = sum(s1)/N;
                gn_chi1(1,n2,n,l+2) =(16/27)*avgF1*Volume*(gamma^2*P0^3);
                gn_chi1(n2,1,n,l+2)=gn_chi1(1,n2,n,l+2);
                
            end
            if( n-2*n2+1-l<0 && n-n2+1-l<=n_ch && 2<=n-n2+1-l && n2-n+l~=0 && 1~=n2)
                w0 = R(1,:)+R(2,:)-R(3,:)+l*2*pi*q;%*l*50e9=*2*pi*q
                arg1 = (R(1,:)-R(3,:)+((1-l-n2)*2*pi*q)).*(R(2,:)-R(3,:)+((n2-n)*2*pi*q));%*50e9=2*pi*q
                argPD1 = arg1;
                ss1 = exp(1i*argPD1*PD).*(exp(1i*beta2*arg1*L-alpha*L)-1)...
                    ./(1i*beta2*arg1-alpha).*(w0<pi).*(w0>-pi);%-16e9=pi
                s1 = abs(ss1.*(1-exp(1i*Nspan*arg1*beta2*L))...
                    ./(1-exp(1i*arg1*beta2*L))).^2/Volume;
                avgF1 = sum(s1)/N;
                gn_chi1(n-n2+1-l,n2,n,l+2) =(16/27)*avgF1*Volume*(gamma^2*P0^3);
                gn_chi1(n2,n-n2+1-l,n,l+2)=gn_chi1(n-n2+1-l,n2,n,l+2);
            end
            
            if(1<=n-n2+1-l+n_ch && n-n2+1-l+n_ch <=n_ch && n-2*n2+l-l+n_ch<0 && n2-n+l~=0 && 1~=n2)
                gn_chi1(n-n2+1-l+n_ch,n2,n,l+2)=0;
                gn_chi1(n2,n-n2+1-l+n_ch,n,l+2)=gn_chi1(n-n2+1-l+n_ch,n2,n,l+2);
            end
        end
        
        
    end
end

% EGN term f1
%     % for n=1:n_ch
for l=-1:1
    for n1=1:n_ch
        for n2=1:n_ch
            if (1<=n1+n2-1+l && n1+n2-1+l<=n_ch && n2-1+l==0 && n1~=n2)
                
                w0 = R(1,:)+R(2,:)-R(3,:)+l*2*pi*q;%*l*50e9=*2*pi*q
                arg1 = (R(1,:)-R(3,:)+((n1-1)*2*pi*q)).*(R(2,:)-R(3,:)+((n2-1)*2*pi*q));%*50e9=2*pi*q
                argPD1 = arg1;
                ss1 = exp(1i*argPD1*PD).*(exp(1i*beta2*arg1*L-alpha*L)-1)...
                    ./(1i*beta2*arg1-alpha).*(w0<pi).*(w0>-pi);%-16e9=pi
                s1 = abs(ss1.*(1-exp(1i*Nspan*arg1*beta2*L))...
                    ./(1-exp(1i*arg1*beta2*L))).^2/Volume;
                avgF1 = sum(s1)/N;
                
                chi11(n1,n2,1,l+2) = (16/27)*avgF1*Volume*(gamma^2*P0^3);%-cof
                chi11(n2,n1,1,l+2)=chi11(n1,n2,1,l+2);
                
                w3p =R(2,:)+R(4,:)-R(3,:)+2*l*pi*q;
                arg2 = (R(2,:)-R(3,:)+((n2-1)*2*pi*q)).*(R(4,:)-R(3,:)+((n1-1)*2*pi*q));%(n2-n)
                argPD2 = arg2;
                ss2 = exp(-1i*argPD2*PD).*(exp(-1i*beta2*arg2*L-alpha*L)-1)...
                    ./(-1i*beta2*arg2-alpha).*(w3p<pi).*(w3p>-pi);
                s2 = (1-exp(1i*Nspan*arg1*beta2*L))./(1-exp(1i*arg1*beta2*L)).*ss1...
                    .*(1-exp(-1i*Nspan*arg2*beta2*L))...
                    ./(1-exp(-1i*arg2*beta2*L)).*ss2/Volume;
                avgF2 = real(sum(s2))/N;
                chi21(n1,n2,1,l+2) = (40/81)*avgF2*Volume*(gamma^2*P0^3);
                chi21(n2,n1,1,l+2)=chi21(n1,n2,1,l+2) ;
            end
            
        end
    end
end
% end
%
% 
for l=-1:1
    for n=2:n_ch
        chi11(2:n_ch,2:n_ch,n,l+2)=chi11(1:n_ch-1,1:n_ch-1,n-1,l+2);
        chi21(2:n_ch,2:n_ch,n,l+2)=chi21(1:n_ch-1,1:n_ch-1,n-1,l+2);
        for n2=1:n_ch
            if (1<=1+n2-n+l && 1+n2-n+l<=n_ch && n2-n+l==0 && 1~=n2)
                w0 = R(1,:)+R(2,:)-R(3,:)+l*2*pi*q;%*l*50e9=*2*pi*q
                arg1 = (R(1,:)-R(3,:)+((1-n)*2*pi*q)).*(R(2,:)-R(3,:)+((n2-n)*2*pi*q));%*50e9=2*pi*q
                argPD1 = arg1;
                ss1 = exp(1i*argPD1*PD).*(exp(1i*beta2*arg1*L-alpha*L)-1)...
                    ./(1i*beta2*arg1-alpha).*(w0<pi).*(w0>-pi);%-16e9=pi
                s1 = abs(ss1.*(1-exp(1i*Nspan*arg1*beta2*L))...
                    ./(1-exp(1i*arg1*beta2*L))).^2/Volume;
                avgF1 = sum(s1)/N;
                chi11(1,n2,n,l+2) = (16/27)*avgF1*Volume*(gamma^2*P0^3);%-cof
                chi11(n2,1,n,l+2)=chi11(1,n2,n,l+2);
                w3p = R(2,:)+R(4,:)-R(3,:)+2*l*pi*q;
                arg2 = (R(2,:)-R(3,:)+((n2-n)*2*pi*q)).*(R(4,:)-R(3,:)+((1-n)*2*pi*q));
                argPD2 = arg2;
                ss2 = exp(-1i*argPD2*PD).*(exp(-1i*beta2*arg2*L-alpha*L)-1)...
                    ./(-1i*beta2*arg2-alpha).*(w3p>-pi).*(w3p<pi);
                s2 = (1-exp(1i*Nspan*arg1*beta2*L))./(1-exp(1i*arg1*beta2*L)).*ss1...
                    .*(1-exp(-1i*Nspan*arg2*beta2*L))...
                    ./(1-exp(-1i*arg2*beta2*L)).*ss2/Volume;
                avgF2 = real(sum(s2))/N;
                chi21(1,n2,n,l+2) = (40/81)*avgF2*Volume*(gamma^2*P0^3);
                chi21(n2,1,n,l+2)=chi21(1,n2,n,l+2);
                
            end
            
        end
        
        for n1=1:n_ch
            if (1<=n1+1-n+l && n1+1-n+l<=n_ch && 1-n+l==0 && n1~=1)
                w0 = R(1,:)+R(2,:)-R(3,:)+l*2*pi*q;%*l*50e9=*2*pi*q
                arg1 = (R(1,:)-R(3,:)+((n1-n)*2*pi*q)).*(R(2,:)-R(3,:)+((1-n)*2*pi*q));%*50e9=2*pi*q
                argPD1 = arg1;
                ss1 = exp(1i*argPD1*PD).*(exp(1i*beta2*arg1*L-alpha*L)-1)...
                    ./(1i*beta2*arg1-alpha).*(w0<pi).*(w0>-pi);%-16e9=pi
                s1 = abs(ss1.*(1-exp(1i*Nspan*arg1*beta2*L))...
                    ./(1-exp(1i*arg1*beta2*L))).^2/Volume;
                avgF1 = sum(s1)/N;
                chi11(n1,1,n,l+2) = (16/27)*avgF1*Volume*(gamma^2*P0^3);%-cof
                chi11(1,n1,n,l+2)=chi11(n1,1,n,l+2);
                w3p = R(2,:)+R(4,:)-R(3,:)+2*l*pi*q;
                arg2 = (R(2,:)-R(3,:)+((1-n)*2*pi*q)).*(R(4,:)-R(3,:)+((n1-n)*2*pi*q));
                argPD2 = arg2;
                ss2 = exp(-1i*argPD2*PD).*(exp(-1i*beta2*arg2*L-alpha*L)-1)...
                    ./(-1i*beta2*arg2-alpha).*(w3p>-pi).*(w3p<pi);
                s2 = (1-exp(1i*Nspan*arg1*beta2*L))./(1-exp(1i*arg1*beta2*L)).*ss1...
                    .*(1-exp(-1i*Nspan*arg2*beta2*L))...
                    ./(1-exp(-1i*arg2*beta2*L)).*ss2/Volume;
                avgF2 = real(sum(s2))/N;
                chi21(n1,1,n,l+2) = (40/81)*avgF2*Volume*(gamma^2*P0^3);
                chi21(1,n1,n,l+2)=chi21(n1,1,n,l+2);
            end
        end
        
    end
    
end

NLIN_var_12=chi11+((kur-2)*chi21);

% EGN term f3
for n=1:n_ch
    for l=-1:1
        for n1=1:n_ch
            if (1<=2*n1-n+l && 2*n1-n+l<=n_ch && 2*n1-n+l~=n1)
                
                w0 = R(1,:)+R(2,:)-R(3,:)+l*2*pi*q;%*l*50e9=*2*pi*q
                arg1 = (R(1,:)-R(3,:)+((n1-n)*2*pi*q)).*(R(2,:)-R(3,:)+((n1-n)*2*pi*q));%*50e9=2*pi*q
                argPD1 = arg1;
                ss1 = exp(1i*argPD1*PD).*(exp(1i*beta2*arg1*L-alpha*L)-1)...
                    ./(1i*beta2*arg1-alpha).*(w0<pi).*(w0>-pi);%-16e9=pi
                s1 = abs(ss1.*(1-exp(1i*Nspan*arg1*beta2*L))...
                    ./(1-exp(1i*arg1*beta2*L))).^2/Volume;
                avgF1 = sum(s1)/N;
                chi13(n1,n1,n,l+2) = (16/27)*avgF1*Volume*(gamma^2*P0^3);
                
                w3p = R(1,:)+R(4,:)-R(3,:)+2*l*pi*q;
                arg2 = (R(1,:)-R(3,:)+((n1-n)*2*pi*q)).*(R(4,:)-R(3,:)+((n1-n)*2*pi*q));
                argPD2 = arg2;
                ss2 = exp(-1i*argPD2*PD).*(exp(-1i*beta2*arg2*L-alpha*L)-1)...
                    ./(-1i*beta2*arg2-alpha).*(w3p>-pi).*(w3p<pi);
                s2 = (1-exp(1i*Nspan*arg1*beta2*L))./(1-exp(1i*arg1*beta2*L)).*ss1...
                    .*(1-exp(-1i*Nspan*arg2*beta2*L))...
                    ./(1-exp(-1i*arg2*beta2*L)).*ss2/Volume;
                avgF2 = real(sum(s2))/N;
                chi23(n1,n1,n,l+2) = (16/81)*avgF2*Volume*(gamma^2*P0^3);
                
            end
        end
    end
end

NLIN_var_3=chi13+((kur-2)*chi23);

% EGN_term SCI&X4
for n=1:n_ch
    for l=-1:1
        for n1=1:n_ch
            if (1<=2*n1-n+l && 2*n1-n+l<=n_ch && 2*n1-n+l==n1)
             
                w0 = R(1,:)+R(2,:)-R(3,:)+l*2*pi*q;%*l*50e9=*2*pi*q
                arg1 = (R(1,:)-R(3,:)+((n1-n)*2*pi*q)).*(R(2,:)-R(3,:)+((n1-n)*2*pi*q));%*50e9=2*pi*q
                argPD1 = arg1;
                ss1 = exp(1i*argPD1*PD).*(exp(1i*beta2*arg1*L-alpha*L)-1)...
                    ./(1i*beta2*arg1-alpha).*(w0<pi).*(w0>-pi);%-16e9=pi
                s1 = abs(ss1.*(1-exp(1i*Nspan*arg1*beta2*L))...
                    ./(1-exp(1i*arg1*beta2*L))).^2/Volume;
                avgF1 = sum(s1)/N;
                chi14(n1,n1,n,l+2) = (16/27)*avgF1*Volume*(gamma^2*P0^3);
                
                w3p = R(1,:)+R(4,:)-R(3,:)+2*l*pi*q;
                arg2 = (R(1,:)-R(3,:)+((n1-n)*2*pi*q)).*(R(4,:)-R(3,:)+((n1-n)*2*pi*q));
                argPD2 = arg2;
                ss2 = exp(-1i*argPD2*PD).*(exp(-1i*beta2*arg2*L-alpha*L)-1)...
                    ./(-1i*beta2*arg2-alpha).*(w3p>-pi).*(w3p<pi);
                s2 = (1-exp(1i*Nspan*arg1*beta2*L))./(1-exp(1i*arg1*beta2*L)).*ss1...
                    .*(1-exp(-1i*Nspan*arg2*beta2*L))...
                    ./(1-exp(-1i*arg2*beta2*L)).*ss2/Volume;
                avgF2 = real(sum(s2))/N;
                
                w33p = R(1,:)+R(2,:)-R(4,:)+2*l*pi*q;
                arg22 = (R(1,:)+R(2,:)-R(4,:)-R(3,:)+((n1-n)*2*pi*q)).*(R(4,:)-R(3,:)+((n1-n)*2*pi*q));
                argPD22 = arg22;
                ss22 = exp(-1i*argPD22*PD).*(exp(-1i*beta2*arg2*L-alpha*L)-1)...
                    ./(-1i*beta2*arg2-alpha).*(w33p>-pi).*(w33p<pi);
                s22 = (1-exp(1i*Nspan*arg1*beta2*L))./(1-exp(1i*arg1*beta2*L)).*ss1...
                    .*(1-exp(-1i*Nspan*arg2*beta2*L))...
                    ./(1-exp(-1i*arg2*beta2*L)).*ss22/Volume;
                avgF22 = real(sum(s22))/N;
                chi24(n1,n1,n,l+2) = (80/81)*avgF2*Volume*(gamma^2*P0^3)+16/81*avgF22*Volume*(gamma^2*P0^3);
                
                
                w3 = R(5,:)+R(4,:)-R(3,:)+l*2*pi*q;
                arg3 = (R(5,:)-R(3,:)+((n1-n)*2*pi*q)).*(R(4,:)-R(3,:)+((n1-n)*2*pi*q));
                argPD3 = arg3;
                ss3 = exp(-1i*argPD3*PD).*(exp(-1i*beta2*arg3*L-alpha*L)-1)...
                    ./(-1i*beta2*arg3-alpha).*(w3<pi).*(w3>-pi);
                s3 = (1-exp(1i*Nspan*arg1*beta2*L))./(1-exp(1i*arg1*beta2*L)).*ss1...
                    .*(1-exp(-1i*Nspan*arg3*beta2*L))./(1-exp(-1i*arg3*beta2*L)).*ss3;
                chi34(n1,n1,n,l+2) = (16/81)*real(sum(s3.*(w0<pi).*(w0>-pi)))*(gamma^2*P0^3)./N;
                
            end
        end
    end
end

NLIN_var_4=chi14+((kur-2)*chi24)+((kur3-9*kur+12)*chi34);
NLIN_var=NLIN_var_12+NLIN_var_3+NLIN_var_4+gn_chi1;

end
