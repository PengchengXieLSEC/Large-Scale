%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn

function [xreturn, fval, temp, k]=main_algorithm(f, x0, verbose)
temp = [];
n = size(x0, 1);
fwp = @(x)(f(x, 1));
f = @(x)(f(x,0));

Delta = 1; % Initial trust-region radius
Delta_low = 1e-3; % Lower bound on trust-region radius
Delta_upper = 1e3; % Upper bound on trust-region radius
gamma1 = 2; % Increase factor for trust-region radius
gamma2 = 1/5; % Decrease factor for trust-region radius
eta = 1/4; % Threshold for very successful step
eta0 = 1/10; % Threshold for successful step

k = 1; % iteration

K=min(100*n,100000); % total iteration
% xs = zeros(n, K);

% Step 0: Initialization Obtain ya, yb, yc
if verbose
    disp(["Iter:", k])
end
d0 = [1; zeros(n-1, 1)];
ya = x0;
yb = x0 + Delta * d0;
fya=fwp(ya);
fyb=fwp(yb);

if fya >= fyb
    yc = ya + 2 * Delta * d0;
else
    yc = ya - Delta * d0;
end

fyc=fwp(yc);

Y = [ya, yb, yc];
Yvalues = [fya, fyb, fyc];
[~, min_index] = min(Yvalues);
ymin1 = Y(:, min_index(1));

if min_index==1
    fymin1=fya;
elseif min_index==2
    fymin1=fyb;
else
    fymin1=fyc;
end
    


% xs(:, k) = ymin1;
xk = ymin1;
if verbose
    disp(["xk", xk'])
end

[~, max_index] = max(Yvalues);
ymax1 = Y(:, max_index(1));

if norm(ymin1 - ymax1)==0
    ymax1 = Y(:, 3-max_index(1)+1);
end

dd1 = (ymin1 - ymax1) / norm(ymin1 - ymax1); 
dd2 = zeros(n, 1);
t = mod(k, n)+1;
dd2(t) = 1;

if (dot(dd1, dd2) ~= 0)
    for index = 1:n
        dd2 = zeros(n, 1);
        t = t + 1;
        dd2(mod(t,n)+1) = 1;
        if (dot(dd1, dd2) == 0)
            break
        end
    end
end

%option for a random d2

%dd2 = null(dd1');
%dd2 = dd2(:, randi(n-1));




alpha_a = dot(ya - ymin1, dd1);
alpha_b = dot(yb - ymin1, dd1);
alpha_c = dot(yc - ymin1, dd1);



A = [alpha_a, alpha_a^2; alpha_b, alpha_b^2; alpha_c, alpha_c^2;];
b = [fya - fymin1; fyb - fymin1; fyc - fymin1];
tempRes = A \ b;
a_value = tempRes(1);
b_value = tempRes(2);

% Loop start
while true
    % Step 1: Construct the interpolation set
    
    
    
    
    if k>1 && norm(xk-y1)==0
        fxk=fy1;
    elseif k>1 && norm(xk-y2)==0
        fxk=fy2;
    elseif k>1 && norm(xk-y3)==0
        fxk=fy3;
    elseif k>1 && norm(xk-xkm1)==0
        fxk=fxkm1;
    elseif exist('y6','var')==1  && norm(xk-y6)==0
        fxk=fy6;
    elseif exist('y8','var')==1  && norm(xk-y8)==0
        fxk=fy8;
    elseif k>1 && norm(xk-xk_plus)==0
        fxk=fxk_plus;
    elseif exist('xk_plus_again','var')==1  && norm(xk-xk_plus_again)==0
        fxk=fxk_plus_again;
    else
        fxk=fwp(xk);
    end
    
   
    
    
    y1 = xk + Delta * dd2;
    
    
    
    if k>1 && norm(y1-xkm1)==0
        fy1=fxkm1;
    else
        fy1=fwp(y1);
    end
    
    
    
    if fy1 <= fxk
        y2 = xk + 2 * Delta * dd2;
    else
        y2 = xk - Delta * dd2;
    end
    

  if k>1 &&  norm(y2-xkm1)==0
        fy2=fxkm1;
  else
    fy2=fwp(y2);
    end
    
    
    
    
    Y2 = [y1, y2];
    Yvalues2 = [fy1, fy2];
    [~, min_index] = min(Yvalues2);
    ymin2 = Y2(:, min_index(1));
    y3 = ymin2 + Delta * dd1;

    
    if k>1 && norm(y3-xkm1)==0
        fy3=fxkm1;
    else
    fy3=fwp(y3);
    end
    
    
    
    flag = true;
    while(true)
%         tp1 = [ya, yb, yc, y1, y2, y3];
%         tp2 = [fya, fyb, fyc, fy1, fy2, fy3];
       
        % Perform the interpolation for Q_k
        
        temp_1=y1 - xk;
        temp_2=y2 - xk;
        temp_3=y3 - xk;
        
        
        alpha1 = dot(temp_1, dd1);
        beta1 = dot(temp_1, dd2);
        
        alpha2 = dot(temp_2, dd1);
        beta2 = dot(temp_2, dd2);
        
        alpha3 = dot(temp_3, dd1);
        beta3 = dot(temp_3, dd2);
        
        A = [beta1, beta1^2, alpha1 * beta1;
            beta2, beta2^2, alpha2 * beta2;
            beta3, beta3^2, alpha3 * beta3];
        b = [fy1 - fxk - a_value * alpha1 - b_value * alpha1 ^ 2;
            fy2 - fxk - a_value * alpha2 - b_value * alpha2 ^ 2;
            fy3 - fxk - a_value * alpha3 - b_value * alpha3 ^ 2];
        tempRes = A \ b;
        c_value = tempRes(1);
        d_value = tempRes(2);
        e_value = tempRes(3);

        
        % Step 3: Trust-region trial step
        g = [a_value; c_value];
        H = [2 * b_value, e_value; e_value, 2 * d_value];

        
        
        cd 'trust_sub'
        [s, ~] = trust_sub(g, H, Delta);
        xk_plus = xk + [dd1, dd2] * s;
        cd '..'
        fxk_plus=fwp(xk_plus);
        
      
       if norm(s)==0
           
            xkp1 = xk;
            fxkp1=fxk;
                    flag = false;
                break
       else

        rho_k = (fxk_plus - fxk) / (1/2*s'*H*s+g'*s);       
        
        if rho_k >= eta0
            xkp1 = xk_plus;
            fxkp1=fxk_plus;
            temp = [temp;k];
            break;
        else
            
            if (k > 1) && (norm(xk-xkm1)~=0)
                
                temp_4_again=xkm1 - xk;
                temp_5_again=xk_plus - xk;
                
                alpha4_again = dot(temp_4_again, dd1);
                beta4_again = dot(temp_4_again, dd2);
                
                alpha5_again = dot(temp_5_again, dd1);
                beta5_again = dot(temp_5_again, dd2);
                
                A = [alpha1, alpha1 ^ 2, beta1, beta1 ^ 2, alpha1 * beta1;
                    alpha2, alpha2 ^ 2, beta2, beta2 ^ 2, alpha2 * beta2;
                    alpha3, alpha3 ^ 2, beta3, beta3 ^ 2, alpha3 * beta3;
                    alpha4_again, alpha4_again ^ 2, beta4_again, beta4_again ^ 2, alpha4_again * beta4_again;
                    alpha5_again, alpha5_again ^ 2, beta5_again, beta5_again ^ 2, alpha5_again * beta5_again];
               
                
                b = [fy1 - fxk; fy2 - fxk; fy3 - fxk;
                    fxkm1 - fxk; fxk_plus - fxk];
                tempRes = A \ b;
                a_value_again = tempRes(1);
                b_value_again = tempRes(2);
                c_value_again = tempRes(3);
                d_value_again = tempRes(4);
                e_value_again = tempRes(5);
    
                % Step 3: Trust-region trial step
                g_again = [a_value_again; c_value_again];
                H_again = [2 * b_value_again, e_value_again; e_value_again, 2 * d_value_again];
                cd 'trust_sub'
                [s_again, ~] = trust_sub(g_again, H_again, Delta);
                xk_plus_again = xk + [dd1, dd2] * s_again;
                cd '..'
                fxk_plus_again=fwp(xk_plus_again);
                if fxk_plus_again < fxk_plus
                    xk_plus=xk_plus_again;
                    fxk_plus=fxk_plus_again;
                end
                rho_k = (fxk_plus - fxk) / ((1/2*s_again'*H_again*s_again+g_again'*s_again));
                
                
                if rho_k >= 0
                    xkp1 = xk_plus;
                    fxkp1=fxk_plus;
                    temp = [temp;k];
                    break;
                else
                    xkp1 = xk;
                    fxkp1=fxk;
                    flag = false;
                    break
                    % end
                end
                
            else   
                xkp1 = xk;
                fxkp1=fxk;
                    flag = false;
                break
            end
            
        
        end
        
       end
    end
    
    oldDelta = Delta;
    
    % Step 4: Update the trust-region radius and the subspace
    if (Delta < Delta_low) || (k > K)
        xreturn = xk;
        fval = fxk;
        return;
    end
    
    if rho_k >= eta
        %0714-new-added
        if Delta<=Delta_upper
            Delta = gamma1 * Delta;
        else
            Delta = 1 * Delta;
        end
        %0714-new-added
    else
        %if norm(s,2)>1e-6
        Delta = gamma2 * Delta;
        %end
    end

   

    
    if flag
        dd1new = (xkp1 - xk) / norm(xkp1 - xk);
    else
        dd1new = dd1;
    end
    
    tp1 = dot(dd1new, dd1);
    tp2 = dot(dd1new, dd2);
    
    dd2star = [-tp2;tp1];
    dd2star = [dd1, dd2] * dd2star;
    
    % dd2 = null([dd1, dd2]');
    % dd2 = dd2(:, randi(n-2));
    
    
    dd2new = zeros(n, 1);
    t = mod(k, n)+1;
    dd2new(t) = 1;
    
    if (dot(dd1, dd2new) ~= 0 || dot(dd2, dd2new) ~= 0)
        for index = 1:n
            dd2new = zeros(n, 1);
            t = t + 1;
            dd2new(mod(t,n)+1) = 1;
            if (dot(dd1, dd2new) == 0 && dot(dd2, dd2new) == 0)
                break
            end
        end
    end

%option for a random d2

%dd2new = null(dd1');
%dd2new = dd2new(:, randi(n-1));






    dd1old = dd1;
    dd2old = dd2;
    
    
    
    dd1 = dd1new;
    dd2 = dd2new;
    
    
    % Perform the interpolation for Q_k^+
    temp_plus_1=y1 - xkp1;
    temp_plus_2=y2 - xkp1;
    temp_plus_3=y3 - xkp1;
    temp_plus_4=xk - xkp1;
    
    alpha1 = dot(temp_plus_1, dd1);
    beta1 = dot(temp_plus_1, dd2star);
    
    alpha2 = dot(temp_plus_2, dd1);
    beta2 = dot(temp_plus_2, dd2star);
    
    alpha3 = dot(temp_plus_3, dd1);
    beta3 = dot(temp_plus_3, dd2star);
    
    alpha4 = dot(temp_plus_4, dd1);
    beta4 = dot(temp_plus_4, dd2star);
    
    if k < 2
        xtemp = x0;
        fxtemp =f(xtemp); 
    else
        xtemp = xkm1;
        fxtemp = fxkm1;
    end
    
    temp_plus_5=xtemp - xkp1;
    
    alpha5 = dot(temp_plus_5, dd1);
    beta5 = dot(temp_plus_5, dd2star);

    y6=xk+oldDelta*dd1old;
    temp_plus_6=y6 - xkp1;   
    alpha6 = dot(temp_plus_6, dd1);
    beta6 = dot(temp_plus_6, dd2star);
   
    y8 = xk + sqrt(2) / 2 * oldDelta * dd2old + sqrt(2) / 2 * oldDelta * dd1old;
    temp_plus_8=y8 - xkp1;
    alpha8 = dot(temp_plus_8, dd1);
    beta8 = dot(temp_plus_8, dd2star);
    % 0714-new-added
    
    
    A1 = [alpha1, alpha1 ^ 2, beta1, beta1 ^ 2, alpha1 * beta1];
    A2 = [alpha2, alpha2 ^ 2, beta2, beta2 ^ 2, alpha2 * beta2];
    A3 = [alpha3, alpha3 ^ 2, beta3, beta3 ^ 2, alpha3 * beta3];
    A4 = [alpha4, alpha4 ^ 2, beta4, beta4 ^ 2, alpha4 * beta4];
    A5 = [alpha5, alpha5 ^ 2, beta5, beta5 ^ 2, alpha5 * beta5];
    A6 = [alpha6, alpha6 ^ 2, beta6, beta6 ^ 2, alpha6 * beta6];
    A8 = [alpha8, alpha8 ^ 2, beta8, beta8 ^ 2, alpha8 * beta8];
        
        
        
    b1 = fy1 - fxkp1;    
    b2 = fy2 - fxkp1;
    b3 = fy3 - fxkp1;
    b4 = fxk - fxkp1;
    b5 = fxtemp - fxkp1;
    
    
    

        
    Yvaluesfordrop = [fy1, fy2];    
    [~, min_indexfordrop] = min(Yvaluesfordrop);
   
    if det([A1; A2; A3; A4; A5])~=0
       tempRes = [A1; A2; A3; A4; A5] \ [b1; b2; b3; b4; b5];
    elseif det([A1; A6; A3; A4; A5])~=0 && min_indexfordrop(1) == 1
             fy6=fwp(y6);
             b6 = fy6 - fxkp1;
              tempRes = [A1; A6; A3; A4; A5] \ [b1; b6; b3; b4; b5];
    elseif det([A6; A2; A3; A4; A5])~=0 
             fy6=fwp(y6);
             b6 = fy6 - fxkp1;
              tempRes = [A6; A2; A3; A4; A5] \ [b6; b2; b3; b4; b5];
    elseif det([A1; A2; A3; A8; A5])~=0 
              fy8=fwp(y8);
              b8 = fy8 - fxkp1;
            tempRes = [A1; A2; A3; A8; A5] \ [b1; b2; b3; b8; b5];
    else
             fy6=fwp(y6);
             b6 = fy6 - fxkp1;
             fy8=fwp(y8);
             b8 = fy8 - fxkp1;
            tempRes = [A1; A2; A3; A8; A6] \ [b1; b2; b3; b8; b6];
    end
                
           
   
    a_value = tempRes(1);
    b_value = tempRes(2);
    c_value = tempRes(3);
    d_value = tempRes(4);
    e_value = tempRes(5);
   
    
    % increase k
     k = k + 1;
     xkm1=xk;
     fxkm1=fxk;
     xk=xkp1;
     
     
     
     
    %     if mod(k,200)==1
%    if mod(k,50)==1
%     k
%    end
    %     end
    
end % end loop

end






