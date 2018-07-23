struct DataMH{vType<:Real}
	h::Vector{vType}
	m::Vector{vType}
end

struct SwParticle{psiType<:Real}
	ψ::psiType
end

function gr!(storage, x, h::Float64, p::MVector)
    storage[1] = -0.5*sin(psi-x[1])+h*sin(x[1])
end

function apply_field(ψ, ϕ, h)
    energy(x) = 0.5*sin(ψ-x[1])^2- h*cos.(x[1])
    gradient(storage, x) = gr!(storage, x, h, p)
    #optimization_result = optimize(energy, [ϕ], Newton())
    optimization_result = optimize(energy,[ϕ], BFGS())
    return Optim.minimizer(optimization_result)[1]
end

function obtain_hysteresis_loop(p::SwParticle)
	hstep = 0.01
    hmax = 2.0
    h = [collect(hmax:-hstep:-hmax); collect(-hmax:hstep:hmax)]
    len = length(h)
    m = zeros(len)
    
    ϕ = p.ψ
    for i=1:len
        ϕ = apply_field(p.ψ, ϕ, h[i])
        m[i] = cos(ϕ)
    end
	return DataMH(h, m)
end

function calculate_particles(psis, h)
    n = length(psis)
    len = length(h)
    phis = copy(psis)
    
    m = zeros(len, n)
    for p = 1:n
        for i = 1:len
            phis[p] = apply_field(psis[p], phis[p], h[i])
            m[i,p] = cos(phis[p])
        end
    end
    return m
end

@recipe function plot_one_sw_particle(data::DataMH)
	legend --> :none
	xlabel --> "Field"
	ylabel --> "Magnetization"
	xlim --> (-2.2, 2.2)
	ylim --> (-1.1,1.1)
    aspect_ratio --> [1 2]
    
    @series begin
        x = vcat((-2.2:0.01:2.2))
        y = vcat(zeros(441))
        linecolor --> :black
        linestyle --> :dash
        (x, y)
    end
    
    @series begin
        x = vcat(zeros(221))
        y = vcat((-1.1:0.01:1.1))
        linecolor --> :black
        linestyle --> :dash
        (x, y)
	end
	
	@series begin
        x = vcat((-1:0.01:1), ones(201), (1:-0.01:-1), -1*ones(201))
        y = vcat(ones(201), (1:-0.01:-1), -1*ones(201), (-1:0.01:1))
        linecolor --> :red
        linestyle --> :dash
        (x, y)
	end
	
    @series begin
        linecolor --> :blue
        linewidth --> 2
        (data.h, data.m)
	end
end
