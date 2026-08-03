# derivative_constants.jl
module DerivativeConstants

export set_constants

function set_constants(dx, dy, dz)
    return (
        co_dx  = 1 / (2 * dx),
        co_dy  = 1 / (2 * dy),
        co_dz  = 1 / (2 * dz),
        co_dxx = 1 / dx^2,
        co_dyy = 1 / dy^2,
        co_dzz = 1 / dz^2,
        co_dxy = 1 / (4.0 * dx * dy),
        co_dxz = 1 / (4.0 * dx * dz),
        co_dyz = 1 / (4.0 * dy * dz)
    )
end

end
