module MagneticParticles

using Optim, StaticArrays, RecipesBase

include("sw_particle.jl")
include("random_sphere_points.jl")
include("usov_particle.jl")


export SwParticle, DataMH
export obtain_hysteresis_loop, calculate_particles

export get_points_spherical_random
export get_points_cartesian_random
export get_points_spherical_evenly
export random_sphere_pack
export create_for_constant_concentration, apply_field!, draw_representation, magnetize_particle

end