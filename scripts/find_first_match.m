function fname = find_first_match(d, p)
% find_first_match  Return first matching file path or empty string
%   fname = find_first_match(directory, pattern)
list = dir(fullfile(d, p));
if ~isempty(list)
    fname = fullfile(d, list(1).name);
else
    fname = '';
end
end