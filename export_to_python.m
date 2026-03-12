% helper function to export the computed u_opt and certificate_opt from the SOP design to a Python file that can be read in python for verification and testing
function export_to_python(u_opt, certificate_opt, params, file_name)
    % u_opt: symbolic expression for the computed controller with bounded control inputs, should be a symbolic vector of size (m x 1) where m is the number of control inputs
    % certificate_opt: symbolic expression for the computed certificate with bounded control inputs, should be a symbolic expression
    % params: struct containing all the parameters for export
    % file_name: name of the file to export the results, should be a string ending with .py

    % Convert symbolic expressions to Python-compatible strings
    u_opt_str = cell(size(u_opt));

    for i = 1:length(u_opt)
        u_opt_str{i} = matlab_expr_to_python(char(u_opt(i)));
    end

    certificate_opt_str = matlab_expr_to_python(char(certificate_opt));

    % Collect all symbolic variables across all expressions to declare them
    all_vars = symvar([u_opt(:); certificate_opt]);
    var_names = arrayfun(@char, all_vars, 'UniformOutput', false);

    % Append timestamp to file name (before the extension)
    [dir, base, ext] = fileparts(file_name);
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    stamped_name = fullfile(dir, sprintf('%s_%s%s', base, timestamp, ext));

    % write the results to a python file
    fid = fopen(stamped_name, 'w');

    % Write sympy import and symbol definitions so expressions are immediately usable
    fprintf(fid, 'from sympy import *\n\n');
    fprintf(fid, '# Symbolic variables\n');
    fprintf(fid, '%s = symbols(''%s'')\n\n', strjoin(var_names, ', '), strjoin(var_names, ' '));

    % Write expressions
    fprintf(fid, 'u_opt = [\n');

    for i = 1:length(u_opt_str)
        fprintf(fid, '    %s,\n', u_opt_str{i});
    end

    % write the certificate expression
    fprintf(fid, ']\n\n');
    fprintf(fid, 'certificate_opt = %s\n', certificate_opt_str);

    % write all the parameters in the params struct to the python file as comments for reference
    fprintf(fid, '\n\n');
    fprintf(fid, '# Parameters\n');

    for field = fieldnames(params)'
        val = params.(field{1});

        if isa(val, 'sym')

            if isscalar(val)
                val_str = matlab_expr_to_python(char(val));
            else
                % symbolic vector/matrix: convert each element, output as nested Python list
                [nr, nc] = size(val);
                rows = cell(nr, 1);

                for r = 1:nr
                    row_strs = arrayfun(@(v) matlab_expr_to_python(char(v)), val(r, :), 'UniformOutput', false);
                    rows{r} = ['[' strjoin(row_strs, ', ') ']'];
                end

                val_str = ['[' strjoin(rows, ', ') ']'];
            end

        else
            % numeric scalar or array: use mat2str for faithful representation
            val_str = mat2str(val);
        end

        fprintf(fid, '# %s = %s\n', field{1}, val_str);
    end

    fclose(fid);
end

function py_str = matlab_expr_to_python(matlab_str)
    % Convert a MATLAB symbolic expression string to a Python/sympy-compatible string
    py_str = matlab_str;
    % MATLAB uses ^ for exponentiation; Python uses **
    py_str = strrep(py_str, '^', '**');
end
