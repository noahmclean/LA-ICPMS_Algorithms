%% Matrix version of Levenberg-Marquardt

function [params, SSE] = LevenbergMarquardt_scalar_v2(fitfun, t, y, p0, maxiter, chiTolerance)

lambda = 1000;

switch fitfun
    case 'exp' 
        a = p0(1); b = p0(2); c = p0(3); m = 3;
        yhat = a*exp(b*t) + c;
end

n = length(t); %m = length(p0);
onesv = ones(n,1);

r = y - yhat;
L = r'*r;
%disp(['  L0 = ' num2str(L)])

G = zeros(m,1);
H = zeros(m,m);
iters = 0;
while iters < maxiter
    
    switch fitfun

        case 'exp'

            Ja =  exp(b*t);
            Jb  = a*t.*Ja;  %scalar multiplication
            Jc =  onesv;
            Jabc = [Ja Jb Jc]; 
            
            G = -Jabc'*r;
            
            H = Jabc'*Jabc;

    end %switch fitcase
    
    h = -(H + lambda*diag(diag(H))) \ G; %calculate update
    
    switch fitfun %for trial update

        case 'exp'
            a_new = a + h(1);
            b_new = b + h(2);
            c_new = c + h(3);
            yhat = a_new*exp(b_new*t) + c_new;
    end
    
    rnew = y - yhat; %test r
    Lnew = rnew'*rnew; %test L
    
    if Lnew > L  %if things got worse, try again
        
        lambda = lambda * 10;
        iters = iters + 1;
        
        %switch fitfun
        %    case 'exp',  disp(['Lnew = ' num2str(Lnew) ' rejected, ' num2str([a b c])])
        %end
        %next up is top of while
        
    else % if things got better

    switch fitfun
            
        case 'exp'
            a = a_new;
            b = b_new;
            c = c_new;
            %disp(['Lnew = ' num2str(Lnew) ' accepted, ' num2str(a) ' ' num2str(b) ' ' num2str(c)])
            
    end %fitcase, update on success

        if (1 - Lnew/L) < chiTolerance  %if not(Lnew>L) and solved
            
            switch fitfun
                case 'exp'
                    params = [a b c];
            end %switch fitcase
            SSE = rnew'*rnew;
            %disp('scalar solved.')  
            return
            
        else %accept function-independent values, lambda /= 10
            
            lambda = lambda/10;
            L = Lnew;
            r = rnew;
            iters = iters + 1;
            
        end %if solved
        
    end %if things are better/worse
    
end % while

%% If failure, try again with -a0 and -b0

lambda = 1000;

switch fitfun
    case 'exp' 
        a = -p0(1); b = -p0(2); c = p0(3); m = 3;
        yhat = a*exp(b*t) + c;
end

n = length(t); %m = length(p0);
onesv = ones(n,1);

r = y - yhat;
L = r'*r;
%disp(['  L0 = ' num2str(L)])

G = zeros(m,1);
H = zeros(m,m);
iters = 0;
while iters < maxiter
    
    switch fitfun

        case 'exp'

            Ja =  exp(b*t);
            Jb  = a*t.*Ja;  %scalar multiplication
            Jc =  onesv;
            Jabc = [Ja Jb Jc]; 
            
            G = -Jabc'*r;
            
            H = Jabc'*Jabc;

    end %switch fitcase
    
    h = -(H + lambda*diag(diag(H))) \ G; %calculate update
    
    switch fitfun %for trial update

        case 'exp'
            a_new = a + h(1);
            b_new = b + h(2);
            c_new = c + h(3);
            yhat = a_new*exp(b_new*t) + c_new;
    end
    
    rnew = y - yhat; %test r
    Lnew = rnew'*rnew; %test L
    
    if Lnew > L  %if things got worse, try again
        
        lambda = lambda * 10;
        iters = iters + 1;
        
        %switch fitfun
        %    case 'exp',  disp(['Lnew = ' num2str(Lnew) ' rejected, ' num2str([a b c])])
        %end
        %next up is top of while
        
    else % if things got better

    switch fitfun
            
        case 'exp'
            a = a_new;
            b = b_new;
            c = c_new;
            %disp(['Lnew = ' num2str(Lnew) ' accepted, ' num2str(a) ' ' num2str(b) ' ' num2str(c)])
            
    end %fitcase, update on success

        if (1 - Lnew/L) < chiTolerance  %if not(Lnew>L) and solved
            
            switch fitfun
                case 'exp'
                    params = [a b c];
            end %switch fitcase
            SSE = rnew'*rnew;
            %disp('scalar solved.')  
            return
            
        else %accept function-independent values, lambda /= 10
            
            lambda = lambda/10;
            L = Lnew;
            r = rnew;
            iters = iters + 1;
            
        end %if solved
        
    end %if things are better/worse
    
end % while

%if you never get there in maxiters
params = -999;
SSE = -999;
disp('LM failure')
