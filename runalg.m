%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn


function frec = runalg(solver, problem, options)

    global fvals fval_history
    global F_data X_data NF_data prob ip is

    rhoend = 1e-5;
    tol = rhoend;
    maxfun = 10000; 
    ftarget = -Inf;
    randomizex0 = 0;
    noiselevel = 0;
    nr = 10;
    

    if (strcmpi('ALL', problem))
        filename = 'problems';
        fileID = fopen(filename);
        C = textscan(fileID, '%s %f');
        tstn = C{2};
        prob = C{1};
        fclose(fileID);
    else
        prob = {problem};
    end
    np = length(prob);

    if (strcmpi('ALL', solver))
        sol = textread('solvers', '%s');
    else
        sol = {solver};
    end
    ns = length(sol);
    temp_lenfavls_old=zeros(ns,1); 

    F_data = zeros(np, ns);
    X_data = cell(np, 1);
    for iip = 1:np
        fvals{iip} = zeros(1, ns * tstn(ip));
    end
    NF_data = zeros(np, ns);

    frec = NaN(np, ns, nr, 200*maxfun);
    fvals = cell(1, ns);
    for iis = 1:ns
        fvals{iis} = [];
    end
    for ip = 1:np
        disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp(strcat(int2str(ip), '. ', prob{ip}, ':'));
        [x0, ~, ~, ~] = setuptest(prob{ip}, tstn(ip));
        maxfun = 100*length(x0); 
        for is = 1:ns
            solv = sol{is};
            if (np * nr <= 20)
                disp('-------------------------------------------------------');
            end
            display(strcat(sol{is}, ':'));
            for ir = 1:nr

                fval_history = [];
             
                if ip==1
                lenfavls = size(fvals{is}, 1);
                else
                    lenfavls = size(fvals{is}, 1)-temp_lenfavls_old(is,1);
                end
              
                if  strcmp(solv, 'twoDMoSub_4test') 
                    Func_wraper = @(x, flg)(evalobjfun(prob{ip}, x, tstn(ip), flg));
                    solverHub = str2func(solv);
                    [x, ~, ~, ~] = solverHub(Func_wraper, x0, [], tstn(ip));
            end
                if (isnan(x))
                    nf= maxfun; 
                    X_data{ip}(is * tstn(ip) - (tstn(ip) - 1):is * tstn(ip)) = x0';

                       F_data(ip, is)= maxfun; 
                     if ip==1

                    lenfavls = maxfun; %20230731
                     else
                    lenfavls = maxfun; %20230731
                     end
                     
                else
                    
                     if ip==1
                    lenfavls = size(fvals{is}, 1);
                     else
                    lenfavls = size(fvals{is}, 1)-temp_lenfavls_old(is,1);
                     end
                     
                    nf = lenfavls;
                    temp_lenfavls_old(is,1)=temp_lenfavls_old(is,1)+lenfavls;

                    X_data{ip}(is * tstn(ip) - (tstn(ip) - 1):is * tstn(ip)) = x';
                    F_data(ip, is) = evalfun(prob{ip}, x, tstn(ip));
                end
                nf = min(nf, maxfun);
                NF_data(ip, is) = nf;

               
            end
        end
        disp('--------------------------------------------------------------------');

        disp(['Solving problem ' num2str(ip) '  Process: ' num2str(ip / np * 100) '%']);

        disp('--------------------------------------------------------------------');
    end

    writecell(X_data, 'X_data.txt');

end

function f = evalobjfun(fun, x, n, flag)
    global fvals is fval_history;
    
    
    [f, ~] = evalfun(fun, x, n);
    
    
        if flag
         fval_history = [fval_history, f];
         fvals{is} = [fvals{is}; f];
        end
    
    
    
end
