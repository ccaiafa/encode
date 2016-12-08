root = pwd;
addpath(genpath(root));

% Compile Mtransp_times_b function
try 
    %cd('mexfiles')
    %mex -largeArrayDims Mtransp_times_b.c Mtransp_times_b_sub.c -output Mtransp_times_b -DNDEBUG
    mex -largeArrayDims CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' OPTIM='-O3 -DNDEBUG -fopenmp' LDOPTIMFLAGS='-O3' COPTIMFLAGS='-O3 -DNDEBUG -fopenmp' INCLUDE='$INCLUDE -I"/N/soft/rhel6/gcc/4.9.2/lib/gcc/x86_64-unknown-linux-gnu/4.9.2/include"' Mtransp_times_b.c Mtransp_times_b_sub.c -output Mtransp_times_b
   
    
    fprintf('Successfully compiled Mtransp_times_b.\n');
    cd(root)
catch lasterr
    cd(root)
    fprintf('Could not compile Mtransp_times_b.');
    fprintf('WARNING: You can still use the ".m" version but it is extremelly slow!!!!!. You should review how to install an appropriate compiler and use it to generate the mex version on you system. See for example: http://www.mathworks.com/support/compilers ');
    rethrow(lasterr);
end

% Compile M_times_w function
try 
    %cd('mexfiles')
    mex -largeArrayDims CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' OPTIM='-O3 -DNDEBUG -fopenmp' LDOPTIMFLAGS='-O3' COPTIMFLAGS='-O3 -DNDEBUG -fopenmp' INCLUDE='$INCLUDE -I"/N/soft/rhel6/gcc/4.9.2/lib/gcc/x86_64-unknown-linux-gnu/4.9.2/include"' M_times_w.c M_times_w_sub.c -output M_times_w
    
    fprintf('Successfully compiled compute_diag.\n');
    cd(root)
catch lasterr
    cd(root)
    fprintf('Could not compile M_times_w.');
    fprintf('WARNING: You can still use the ".m" version but it is extremelly slow!!!!!. You should review how to install an appropriate compiler and use it to generate the mex version on you system. See for example: http://www.mathworks.com/support/compilers ');
    rethrow(lasterr);
end

% Compile compute_diag function
try 
    %cd('mexfiles')
    mex -largeArrayDims compute_diag.c compute_diag_sub.c -output compute_diag -DNDEBUG
    
    fprintf('Successfully compiled compute_diag.\n');
    cd(root)
catch lasterr
    cd(root)
    fprintf('Could not compile compute_diag.');
    fprintf('WARNING: You can still use the ".m" version but it is extremelly slow!!!!!. You should review how to install an appropriate compiler and use it to generate the mex version on you system. See for example: http://www.mathworks.com/support/compilers ');
    rethrow(lasterr);
end


