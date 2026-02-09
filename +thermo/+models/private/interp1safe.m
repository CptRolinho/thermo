function yq = interp1safe(x, v, xq, method, extrap)
%interp1safe Safe 1-D interpolation with input cleaning and bounds checks.
%
% yq = interp1safe(x, v, xq)
% yq = interp1safe(x, v, xq, method)
% yq = interp1safe(x, v, xq, method, extrap)
%
% - Filters out non-finite (NaN/Inf) pairs in (x,v)
% - Requires xq to be finite scalar
% - Sorts and removes duplicate x entries (keeps first occurrence)
% - Errors on out-of-range unless extrap is provided
% NOTE:
%   Any non-finite (NaN/Inf) pairs in the tabulated data (x,v) are filtered
%   out before interpolation. This makes the function robust to incomplete
%   textbook tables while keeping strict bounds checking for xq.

    if nargin < 4 || isempty(method), method = 'linear'; end
    if nargin < 5, extrap = []; end

    validateattributes(xq, {'numeric'}, {'scalar','finite'}, mfilename, 'xq');

    % Force vectors
    x = x(:);
    v = v(:);

    if numel(x) ~= numel(v)
        error('interp1safe:SizeMismatch', 'x and v must have the same number of elements.');
    end

    % Filter non-finite pairs (common in tables with missing entries)
    m = isfinite(x) & isfinite(v);
    x = x(m);
    v = v(m);

    if numel(x) < 2
        error('interp1safe:InsufficientData', ...
            'Not enough finite data points after filtering (need >=2).');
    end

    % Sort by x
    [x, idx] = sort(x);
    v = v(idx);

    % Remove duplicate x values (interp1 requires monotonic/unique grid)
    [xU, iu] = unique(x, 'stable');
    vU = v(iu);

    % Range check (unless extrap specified)
    xmin = xU(1); xmax = xU(end);
    if isempty(extrap)
        if xq < xmin || xq > xmax
            error('interp1safe:OutOfRange', ...
                'Query xq=%.8g outside data range [%.8g, %.8g].', xq, xmin, xmax);
        end
        yq = interp1(xU, vU, xq, method);
    else
        yq = interp1(xU, vU, xq, method, extrap);
    end

    if ~isfinite(yq)
        error('interp1safe:NonFiniteResult', ...
            'Interpolation returned non-finite result at xq=%.8g.', xq);
    end
end
