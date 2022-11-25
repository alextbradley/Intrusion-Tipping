function dhdx = compute_slope(h,dx)
% compute the slope dh/dx of h using second order centreded FD in the interior and one sided FD at edges. 
dhdx = zeros(1, length(h));
dhdx(2:end-1) = (h(1:end-2) - h(3:end))/2/dx;
dhdx(1)       = -(-3/2 * h(1) + 2*h(2) - 1/2 *h(3))/dx;
dhdx(end)     = -(1/2 * h(end-2) - 2*h(end-1) + 3/2*h(end))/dx;
%dhdx = smooth(dhdx);
end
