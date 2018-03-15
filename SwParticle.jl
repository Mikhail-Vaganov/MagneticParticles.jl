module StonerWohlfarthModel 

using Optim, StaticArrays, Plots


function gr!(storage, x, h::Float64, p::MVector)
    storage[1] = -0.5*sin(psi-x[1])+h*sin(x[1]);
end

function apply_field!(ψ, ϕ, h)
    energy(x) = 0.5*sin(ψ-x[1])^2- h*cos.(x[1]);
    gradient(storage, x) = gr!(storage, x, h, p);
    #optimization_result = optimize(energy, [ϕ], Newton());
    optimization_result = optimize(energy,[ϕ], BFGS());
    return Optim.minimizer(optimization_result)[1];
end

function draw_representation(ψ::Float64)
    hstep = 0.01
    hmax = 2.0
    h = [collect(hmax:-hstep:-hmax); collect(-hmax:hstep:hmax)]
    len = length(h)
    println(len)
    m = zeros(len)
    
    ϕ = ψ
    for i=1:len
        ϕ = apply_field!(ψ, ϕ, h[i])
        m[i] = cos(ϕ)
    end
    plot(h, m, xlim = (-2.2, 2.2), ylim = (-1.1,1.1), aspect_ratio=[1 2], legend=:none, xlabel="Field", ylabel = "Magnetization");
end

function calculate_particles(psis::Array{Float64}, h::Array{Float64})
    n = length(psis)
    len = length(h)
    phis = copy(psis)
    
    m = zeros(n, len)
    for p = 1:n
        for i = 1:len
            phis[p] = apply_field!(psis[p], phis[p], h[i])
            m[p,i] = cos(phis[p])
        end
    end
    return m
end

export draw_representation, calculate_particles
end