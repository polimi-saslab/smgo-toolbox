function f = real2normd( x, bnds )
    x_len = size(x, 2);
    f = (x - repmat(bnds(:, 1), 1, x_len)) ./ repmat(bnds(:, 2) - bnds(:, 1), 1, x_len);
end