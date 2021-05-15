"""
linesearch!(x, fvec, g, p, stpmax, func, tolx, args...) -> f
% Computes the optimal step by minimizing the residual sum of squares
%
% INPUTS
%   xold:     actual point
%   fold:     residual sum of squares at the point xold
%   g:        gradient
%   p:        Newton direction
%   stpmax:   maximum step
%   func:     name of the function
%   j1:       equations index to be solved
%   j2:       unknowns index
%   tolx:     tolerance parameter
%   varargin: list of arguments following j2
%
% OUTPUTS
%   x:        chosen point
%   fvec:     residuals vector
%   f:        residual sum of squares value for a given x
%   check=1:  problem of the looping which continues indefinitely
%
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2018 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
"""
function linesearch!(x, fvec, g, p, stpmax, func, tolx, args...)

    alf = 1e-4 
    alam = 1

    nn = length(x)
    summ = sqrt(p'*p)

    if !isfinite(summ)
        error("""
    Some element of Newton direction isn''t finite. Jacobian maybe
    singular or there is a problem with initial values
    """)
    end

    if summ > stpmax
        p = p*stpmax/summ 
    end

    slope = dot(g, p) 

    test = get_test(p, x)
    alamin = tolx/test 

    if alamin > 0.1
        alamin = 0.1
    end

    while 1
        if alam < alamin
            check = 1 
            return (check, f)
        end
        x .+= alam*p
        try
            fvec = func(x,args)
        catch
            fvec = NaN
        end
        f = 0.5*dot(fvec, fvec) 
        if any(isnan.(fvec))
            alam = alam/2 
            alam2 = alam 
            f2 = f 
            fold2 = fold 
        else
            if f <= fold + alf*alam*slope
                check = 0
                break
            else
                if alam == 1
                    tmplam = -slope/(2*(f - fold - slope)) 
                else
                    rhs1 = f - fold - alam*slope 
                    rhs2 = f2 - fold2 - alam2*slope
                    alamsq = alam*alam
                    alam2sq = alam2*alam2
                    a = (rhs1/alamsq - rhs2/alam2sq)/(alam - alam2) 
                    b = (-alam2*rhs1/alamsq + alam*rhs2/alam2sq)/(alam - alam2) 
                    if a == 0
                        tmplam = -slope/(2*b) 
                    else
                        disc = (b^2)-3*a*slope 

                        if disc < 0
                            error("Roundoff problem in linesearch!") 
                        else
                            tmplam = (-b + sqrt(disc))/(3*a) 
                        end
                    end
                    if tmplam > 0.5*alam
                        tmplam = 0.5*alam
                    end
                end
                alam2 = alam 
                f2 = f 
                fold2 = fold 
                alam = max(tmplam, 0.1*alam) 
            end
        end
    end
    check = 0
    return(check, f)
end

function get_test(p, x)
    test = 0.0
    for i in length(p)
        y = abs(p[i])/max(abs(x[i]), 1.0)
        if y > test
            test = y
        end
    end
    return test
end
    
