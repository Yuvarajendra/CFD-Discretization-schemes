function [x, residual, numItr] = gauss_seidel(A, b, x0,ConvCrit)

    % Additional inputs:
   
    rlx=0.1; % Relaxation factor
    n=length(x0);
    x = x0; xold = x0;
    numItr=0;
    residual=1;

    while residual>= ConvCrit % Continue until the residual is smaler than the convergence criteria

        for i=1:n % Loop over gridpoints

            I_lhs=[1:i-1]; % left-hand-side values from the diagonal in each row in A matrix
            I_rhs=[i+1:n]; % right-hand-side values from the diagonal in each row in A matrix
               
            % Compute solution (x) for each grid point (including relaxation): 
            
            x(i) =(1-rlx)*xold(i)+ rlx*(-A(i,I_lhs)*x(I_lhs)/A(i,i) - A(i,I_rhs)*xold(I_rhs)/A(i,i) +  b(i)/A(i,i));
        end
        
        r= b-A*x; % Equation residual

        residual=max(abs(r)); % Maximum value as residual criteria
        numItr=numItr+1; % Number of solver iterations
        xold=x; % Replace old solution with new
    end
end

