function methods = auto_discover_method_dirs(root_dir, required_patterns)
% AUTO_DISCOVER_METHOD_DIRS Recursively discover method data folders under a root path.
% methods(i) fields:
%   - data_dir
%   - method_name
%   - session_name
%   - latest_time

    if nargin < 2 || isempty(required_patterns)
        required_patterns = {'*_PhR_stats.mat'};
    end

    if ~exist(root_dir, 'dir')
        error('Root directory does not exist: %s', root_dir);
    end

    probe_pattern = required_patterns{1};
    probe_files = dir(fullfile(root_dir, '**', probe_pattern));

    if isempty(probe_files)
        methods = struct('data_dir', {}, 'method_name', {}, 'session_name', {}, 'latest_time', {});
        return;
    end

    candidate_dirs = unique(string({probe_files.folder}));
    records = struct('data_dir', {}, 'method_name', {}, 'session_name', {}, 'latest_time', {});

    for i = 1:numel(candidate_dirs)
        folder_i = char(candidate_dirs(i));
        if ~local_has_required_files(folder_i, required_patterns)
            continue;
        end

        [method_name, session_name] = local_extract_names(root_dir, folder_i);

        probe_here = dir(fullfile(folder_i, probe_pattern));
        if isempty(probe_here)
            latest_time = 0;
        else
            latest_time = max([probe_here.datenum]);
        end

        rec.data_dir = folder_i;
        rec.method_name = method_name;
        rec.session_name = session_name;
        rec.latest_time = latest_time;
        records(end+1) = rec; %#ok<AGROW>
    end

    if isempty(records)
        methods = struct('data_dir', {}, 'method_name', {}, 'session_name', {}, 'latest_time', {});
        return;
    end

    all_methods = unique(string({records.method_name}));
    methods = struct('data_dir', {}, 'method_name', {}, 'session_name', {}, 'latest_time', {});

    for k = 1:numel(all_methods)
        m = char(all_methods(k));
        idx = find(strcmp({records.method_name}, m));
        [~, best_local] = max([records(idx).latest_time]);
        methods(end+1) = records(idx(best_local)); %#ok<AGROW>
    end

    [~, order] = sort(lower(string({methods.method_name})));
    methods = methods(order);
end

function tf = local_has_required_files(dir_path, required_patterns)
    tf = true;
    for i = 1:numel(required_patterns)
        if isempty(dir(fullfile(dir_path, required_patterns{i})))
            tf = false;
            return;
        end
    end
end

function [method_name, session_name] = local_extract_names(root_dir, data_dir)
    root_clean = char(string(root_dir));
    data_clean = char(string(data_dir));

    if endsWith(root_clean, filesep)
        root_clean = root_clean(1:end-1);
    end

    rel = strrep(data_clean, [root_clean filesep], '');
    parts = strsplit(rel, filesep);

    if isempty(parts)
        method_name = 'Method';
        session_name = '';
        return;
    end

    session_name = parts{end};

    if numel(parts) >= 2
        method_name = parts{end-1};
    else
        [~, method_name] = fileparts(data_clean);
    end

    if isempty(method_name)
        method_name = 'Method';
    end
end
