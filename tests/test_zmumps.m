function test_zmumps()
    clc;
    disp('==============================================');
    disp('        MUMPS OpenMP Complex Test (ZMUMPS)    ');
    disp('==============================================');
    
    % 1. Check if the MEX file exists
    % Note: When compiling the complex version, the output file is named zmumpsmex.mexw64
    if exist('zmumpsmex', 'file') ~= 3
        error('Error: zmumpsmex.mexw64 file not found! Please run compile_mumps(''z'') first.');
    end
    disp('MEX file (zmumpsmex) check passed.');
    
    % 2. Set the number of OpenMP threads
    n_threads = 4; 
    setenv('OMP_NUM_THREADS', num2str(n_threads));
    fprintf('OMP_NUM_THREADS set to = %d\n', n_threads);
    
    % 3. Prepare complex test data
    n = 5000;
    density = 0.01;
    fprintf('Generating complex test matrix (n=%d, density=%.2f)...\n', n, density);
    
    % Construct complex sparse matrix A = Real + i*Imag
    % Maintain diagonal dominance to ensure numerical stability
    A = (sprand(n, n, density) + speye(n)*10) + ...
        1i * (sprand(n, n, density) + speye(n)*10);
        
    % Right-hand side is also complex
    b = rand(n, 1) + 1i * rand(n, 1);
    
    % 4. Initialize MUMPS
    disp('Initializing ZMUMPS...');
    id = initmumps;
    id.JOB = -1;
    id.SYM = 0; % 0=Unsymmetric (Complex matrices usually use Unsymmetric or Hermitian; using Unsymmetric here)
    
    % [NOTE] The zmumps interface needs to be called here
    % If the zmumps.m wrapper file is not in your folder, you can call the mex file directly:
    if exist('zmumps', 'file') == 2
        id = zmumps(id);
    else
        warning('zmumps.m wrapper file not found, attempting to call zmumpsmex directly.');
        id = zmumpsmex(id);
    end
    
    % 5. Set control parameters
    id.ICNTL(1:4) = 0; % Turn off output
    
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
    
    % 6. Execute the solve
    disp('Starting solve (Complex)...');
    id.JOB = 6; 
    id.RHS = b;
    
    tic;
    if exist('zmumps', 'file') == 2
        id = zmumps(id, A);
    else
        id = zmumpsmex(id, A);
    end
    t_solve = toc;
    
    fprintf('Solve complete! Time taken: %.4f seconds\n', t_solve);
    
    % 7. Verify results
    if ~isempty(id.SOL)
        x = id.SOL;
        % Calculate complex residual
        rel_res = norm(A*x - b) / norm(b);
        fprintf('Relative residual (Complex): %.2e\n', rel_res);
        
        if rel_res < 1e-10
            disp('>> Test successful! Result is correct. <<');
        else
            disp('>> Warning: Residual is too large, result may be incorrect. <<');
        end
    else
        disp('Error: Solution vector not returned.');
    end
    
    % 8. Release memory
    disp('Cleaning up memory...');
    id.JOB = -2;
    if exist('zmumps', 'file') == 2
        id = zmumps(id);
    else
        id = zmumpsmex(id);
    end
    disp('Done.');
end