function xabs = unwrap(x)
% UNWRAP maps position onto absolute coordinates

thresh = 24;        %detect ms jumps by 24 pixels or more (these are typically wrapping errors)
xabs = x;
dx = diff(xabs);
ii = find(abs(dx)>thresh); 
for i=1:numel(ii)
    xabs( (ii(i)+1):end ) = xabs( (ii(i)+1):end ) - dx(ii(i));
end

end