using MagneticParticles

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@test 1==1
@test 1==0
