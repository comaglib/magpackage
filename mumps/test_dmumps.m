function test_dmumps()
    clc;
    disp('==============================================');
    disp('        MUMPS OpenMP Interface Test           ');
    disp('==============================================');
    
    % 1. Check if the MEX file exists
    if exist('dmumpsmex', 'file') ~= 3
        error('Error: dmumpsmex.mexw64 file not found! Please ensure successful compilation and that it is in the current path.');
    end
    disp('MEX file check passed.');
    
    % 2. Set the number of OpenMP threads (Core test point)
    % We set different thread counts to observe performance variation
    n_threads = 4; 
    setenv('OMP_NUM_THREADS', num2str(n_threads));
    fprintf('OMP_NUM_THREADS set to = %d\n', n_threads);
    
    % 3. Prepare test data (Generate a slightly large sparse matrix)
    n = 5000;
    density = 0.01;
    fprintf('Generating test matrix (n=%d, density=%.2f)...\n', n, density);
    A = sprand(n, n, density) + speye(n)*10; % Diagonal dominance to ensure invertibility
    A = A + A'; % Symmetrize (though MUMPS can also handle asymmetric matrices)
    b = rand(n, 1);
    
    % 4. Initialize MUMPS
    disp('Initializing MUMPS...');
    id = initmumps;
    id.JOB = -1;
    id.SYM = 0; % 0=Unsymmetric, 1=Positive Definite Symmetric, 2=General Symmetric (Using 0 for general testing here)
    id = dmumps(id); 
    
    % 5. Set control parameters
    % ICNTL(1)-ICNTL(4) control output verbosity, set to 0 to keep quiet
    id.ICNTL(1:4) = 0; 
    
    % ICNTL(16) is the OpenMP thread count setting, usually set to 0 to follow the OMP_NUM_THREADS environment variable
    % But some versions may require explicit setting
    id.ICNTL(16) = n_threads; 
    
    % Enable output for observation
    id.ICNTL(1) = 6; 
    id.ICNTL(2) = 0; 
    id.ICNTL(3) = 6; 
    id.ICNTL(4) = 2; % Print detailed information
    
    % Set the ordering algorithm
    % 0: AMD, 2: AMF, 3: SCOTCH, 4: PORD, 5: METIS, 7: Auto
    id.ICNTL(7) = 5; % Force use of METIS
    
    % 6. Execute the solve (Analysis + Factorization + Solve)
    disp('Starting solve...');
    id.JOB = 6; 
    id.RHS = b;
    
    tic;
    id = dmumps(id, A);
    t_solve = toc;
    
    fprintf('Solve complete! Time taken: %.4f seconds\n', t_solve);
    
    % 7. Verify results
    if ~isempty(id.SOL)
        x = id.SOL;
        rel_res = norm(A*x - b) / norm(b);
        fprintf('Relative residual: %.2e\n', rel_res);
        
        if rel_res < 1e-10
            disp('>> Test successful! Result is correct. <<');
        else
            disp('>> Warning: Residual is too large, result may be incorrect. <<');
        end
    else
        disp('Error: Solution vector not returned. Check INFO(1) error code.');
        disp(['INFO(1) = ', num2str(id.INFO(1))]);
        disp(['INFO(2) = ', num2str(id.INFO(2))]);
    end
    
    % 8. Release memory
    disp('Cleaning up memory...');
    id.JOB = -2;
    id = dmumps(id);
    disp('Done.');
end