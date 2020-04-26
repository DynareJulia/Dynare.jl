nearbyint(x::Float64) = (abs((x)-floor(x)) < abs((x)-ceil(x)) ? floor(x) : ceil(x))

function get_power_deriv(x::Float64, p::Float64, k::Int64)
  if ( abs(x) < 1e-12 && p > 0 && k > p && abs(p-nearbyint(p)) < 1e-12 )
    return 0.0
  else
      dxp = x^(p-k)
      for i = 1:k
          dxp *= p
          p -= 1
      end
      return dxp
  end
end

