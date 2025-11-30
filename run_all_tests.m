function run_all_tests()
% RUN_ALL_TESTS Automated Test Suite for MagPackage (v2.3 - Stable)
% 
% Features:
%   1. [Progress Bar] Visual waitbar showing test progress.
%   2. [Robustness] Protected waitbar from 'close all' to prevent interruptions.
%   3. [Cleanup] Ensures output buffers and graphical resources are cleared.
%   4. [WARN Detection] Captures and reports warnings from tests.

    % 1. Environment Initialization
    clc;
    fprintf('============================================================\n');
    fprintf('                 MagPackage Automated Test Suite            \n');
    fprintf('============================================================\n');
    
    baseDir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(baseDir, 'src')));
    testDir = fullfile(baseDir, 'tests');
    addpath(testDir);
    
    % Check for MEX files
    if exist('assemble_curl_curl_kernel_mex', 'file') ~= 3
        fprintf('[WARN] MEX files not detected. Attempting to compile...\n');
        try
            make();
        catch
            fprintf(2, '[ERROR] Compilation failed. Aborting tests.\n');
            return;
        end
    end

    % 2. Scan Test Files
    files = dir(fullfile(testDir, 'test_*.m'));
    if isempty(files)
        fprintf('[WARN] No test files found in %s\n', testDir);
        return;
    end
    
    numTests = length(files);
    results = struct('name', {}, 'status', {}, 'time', {}, 'message', {});
    
    fprintf('[INFO] Found %d test scripts. Starting execution...\n', numTests);
    fprintf('============================================================\n');
    
    % --- Initialize Progress Bar ---
    % Set a specific Tag to protect it from 'close all'
    hBar = waitbar(0, 'Initializing test suite...', ...
        'Name', 'MagPackage Test Runner', ...
        'Tag', 'MagPackageWaitbar', ... 
        'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
    setappdata(hBar, 'canceling', 0);
    
    % Register cleanup to ensure bar closes on error/Ctrl+C
    cleanupObj = onCleanup(@() cleanup_resources(hBar));
    
    totalStart = tic;
    
    % 3. Main Test Loop
    for i = 1:numTests
        % Check for user cancellation
        if ~isvalid(hBar) || getappdata(hBar, 'canceling')
            fprintf(2, '\n[USER ABORT] Test suite cancelled by user.\n');
            break;
        end
        
        fileName = files(i).name;
        [~, scriptName] = fileparts(fileName);
        
        % Update Progress Bar
        % Escape underscores for display
        waitbar((i-1)/numTests, hBar, sprintf('Running [%d/%d]: %s', i, numTests, strrep(fileName, '_', '\_')));
        
        fprintf('--> [%02d/%02d] Testing %-35s \n', i, numTests, fileName);
        fprintf('------------------------------------------------------------\n');
        
        % Pre-test Cleanup
        safe_close_all(); % Close other figures, keep waitbar
        lastwarn('');     % Reset warning state
        
        tStart = tic;
        
        % [Core] Run Script with Output Capture
        [txtOutput, success, ME] = capture_execution(scriptName);
        timeSec = toc(tStart);
        
        % Print captured output with indentation
        if ~isempty(txtOutput)
            fprintf('%s', txtOutput);
            if ~endsWith(txtOutput, newline), fprintf('\n'); end
        end
        
        % 4. Determine Status (PASS / FAIL / WARN)
        statusStr = 'PASS';
        msg = '';
        
        if ~success
            statusStr = 'FAIL';
            msg = ME.message;
        else
            [sysWarnMsg, ~] = lastwarn;
            hasCustomWarn = contains(txtOutput, '[WARN]', 'IgnoreCase', false);
            hasSysWarn = ~isempty(sysWarnMsg);
            
            if hasCustomWarn || hasSysWarn
                statusStr = 'WARN';
                if hasSysWarn, msg = sysWarnMsg; 
                else, msg = 'Custom Warning detected.'; end
            end
        end
        
        % 5. Record Results
        results(i).name = fileName;
        results(i).status = statusStr;
        results(i).time = timeSec;
        results(i).message = msg;
        
        % Print Footer
        if strcmp(statusStr, 'FAIL')
            fprintf(2, '>>> RESULT: FAIL (%.2fs) <<<\n', timeSec);
            fprintf(2, '    Error: %s\n', msg);
        elseif strcmp(statusStr, 'WARN')
            fprintf('>>> RESULT: WARN (%.2fs) <<<\n', timeSec);
        else
            fprintf('>>> RESULT: PASS (%.2fs) <<<\n', timeSec);
        end
        fprintf('============================================================\n\n');
        
        drawnow limitrate;
    end
    
    % Finalize Progress Bar
    if isvalid(hBar)
        waitbar(1, hBar, 'Tests Completed. Generating Report...');
        pause(0.5);
    end
    
    totalTime = toc(totalStart);
    
    % Explicitly delete waitbar (onCleanup handles it, but explicit is cleaner here)
    delete(hBar); 
    
    % 6. Generate Summary Report
    print_summary(results, totalTime);
end

function [txt, success, ME] = capture_execution(scriptName)
    % Runs the script in an isolated environment and captures stdout
    success = false;
    ME = [];
    txt = '';
    
    try
        % Use evalc to capture output from the isolated runner
        cmd = sprintf('run_script_isolated(''%s'');', scriptName);
        txt = evalc(cmd);
        success = true;
    catch innerME
        success = false; 
        ME = innerME;
        if isempty(txt)
            txt = '[Output captured stopped due to Error]';
        end
    end
end

function run_script_isolated(scriptName)
    % Wrapper to run script in its own workspace
    eval(scriptName);
end

function safe_close_all()
    % Closes all figures EXCEPT the progress bar
    % This prevents the 'stuck' issue where the loop destroys its own UI
    figs = findall(0, 'Type', 'figure', '-not', 'Tag', 'MagPackageWaitbar');
    if ~isempty(figs)
        close(figs);
    end
end

function cleanup_resources(hBar)
    % Callback for onCleanup
    if isvalid(hBar), delete(hBar); end
end

function print_summary(results, totalTime)
    % Print summary with extra spacing
    fprintf('\n\n'); 
    fprintf('############################################################\n');
    fprintf('                     TEST SUMMARY                           \n');
    fprintf('############################################################\n');
    
    if isempty(results)
        return;
    end
    
    numTests = length(results);
    cntPass = sum(strcmp({results.status}, 'PASS'));
    cntWarn = sum(strcmp({results.status}, 'WARN'));
    cntFail = sum(strcmp({results.status}, 'FAIL'));
    
    fprintf('Total Time: %.2f seconds\n', totalTime);
    fprintf('Total Tests: %d\n', numTests);
    fprintf('  [PASS]: %d\n', cntPass);
    fprintf('  [WARN]: %d\n', cntWarn);
    fprintf('  [FAIL]: %d\n', cntFail);
    fprintf('------------------------------------------------------------\n');
    
    if cntFail > 0 || cntWarn > 0
        fprintf('Attention Required:\n');
        for k = 1:numTests
            st = results(k).status;
            if strcmp(st, 'FAIL')
                fprintf(2, '  [FAIL] %-30s : %s\n', results(k).name, results(k).message);
            elseif strcmp(st, 'WARN')
                fprintf('  [WARN] %-30s : %s\n', results(k).name, results(k).message);
            end
        end
        fprintf('\n[RESULT] SUITE FINISHED WITH ISSUES.\n');
    else
        fprintf('[RESULT] ALL TESTS PASSED PERFECTLY.\n');
    end
    fprintf('############################################################\n');
end