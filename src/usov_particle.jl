
global μ₀=4pi*1e-7

mutable struct UsovParticle
    n::Int
    νs::Array{Float64,2}
    rs2ds::Array{Float64,2}
    ns::Array{Float64,3}
    ψs::Array{Float64,1}
    ϕs::Array{Float64,1}
    λ::Float64
    hprev::Float64
    hsat::Float64

    function UsovParticle(r::Float64,R::Float64,n::Int,η::Float64,K_M2::Float64)
        ps = random_sphere_pack(r, R, n)
        νs = get_points_cartesian_random(n)
        ϕs = acos.(νs[:,3])
        ψs = atan2.(νs[:,2], νs[:,1])
        rs2ds = zeros(n, n)
        ns = zeros(3, n, n)
        for i = 1:n
            for j = 1:n
                ds = ps[j,:] -  ps[i,:]
                rs2ds[i,j] = (r/norm(ds))^3
                ns[:,i,j] = ds/norm(ds)
            end
        end
        λ = μ₀/(6*K_M2)
        return new(n, νs, rs2ds, ns, ψs, ϕs, λ, 0.0, 4.0)
    end
end

function total_energy(n, νs, ns, λ, ϕs, ψs, rs2ds, h)
    es = [ sin.(ϕs).*cos.(ψs) sin.(ϕs).*sin.(ψs)  cos.(ϕs)]
    
    #sum(es.*nus,2) - scallar multiplication of all particle vectors es and nus 
    energy_anis_and_zeeman = sum( 0.5(1-sum(es.*νs, 2).^2) - h*es[:,3]) #external sum - by particles
    energy_dip_dip = 0
    for i = 1:n
        # external sum by j-th neighbouring particles
        #v1 = λ*rs2ds[i+1:n,i]
        #v21 = (es[i,:]'*ns[:,i+1:n,i])
        #v22 = sum( es[i+1:n,:].*ns[:,i+1:n,i]',2 )
        #v2 = v21'.*v22
        #v3 = es[i+1:n,:]*es[i,:]
        # energy_dip_dip -= sum(v1.*(3v2.-v3))
        energy_dip_dip -= sum(λ*rs2ds[i+1:n,i].* (3*(es[i,:]'*ns[:,i+1:n,i])'.*sum( es[i+1:n,:].*ns[:,i+1:n,i]',2 ) - es[i+1:n,:]*es[i,:]) )
    end
    
    return energy_anis_and_zeeman + energy_dip_dip
end

function total_energy_for_loop(n, νs, ns, λ, ϕs, ψs, rs2ds, h)
    es = @. [ sin(ϕs)*cos(ψs) sin(ϕs)*sin(ψs)  cos(ϕs)]
     
    anisotropy = 0.0
    zeeman = 0.0
    dip_dip = 0.0
    for i = 1:n
        anisotropy += 0.5(1-sum(es[i,:].*νs[i,:]).^2)
        zeeman += h*es[i,3]
        for j = i+1:n
            dip_dip += λ*rs2ds[j,i]* (3*sum(es[i,:].*ns[:,j,i])*sum(es[j,:].*ns[:,j,i]) - sum(es[j,:].*es[i,:]))
        end
    end    
    return anisotropy - zeeman - dip_dip
end


function apply_field_by_step!(p::UsovParticle, h::Float64)
    energy(x) = total_energy(p.n, p.νs, p.ns, p.λ, x[1:p.n], x[p.n+1:2p.n], p.rs2ds, h);
    optimization_result = optimize(energy, [p.ϕs; p.ψs], BFGS());
    angles = Optim.minimizer(optimization_result);
    p.ϕs[1:end] = angles[1:p.n]
    p.ψs[1:end] = angles[p.n+1:2p.n]
    p.hprev = h
end

function apply_field!(p::UsovParticle, applied_field::Float64)
    hstep = 0.01
    if applied_field>p.hprev
        field_values = p.hprev:hstep:applied_field
    else
        field_values = p.hprev:-hstep:applied_field
    end
    
    if field_values[end] != applied_field
        push!(field_values,applied_field)
    end
    
    for h in field_values
        apply_field_by_step!(p, h)
    end
        
end

function create_for_constant_concentration(r::Float64,n::Int,η::Float64,K::Float64,M::Float64)
    R = r*(n/η)^(1/3)
    K_M2 = K/(M*M)
    p = UsovParticle(r, R, n, η, K_M2)
end

function draw_representation(p::UsovParticle)
    hstep = 0.01
    hmax = p.hsat
    h = [collect(hmax:-hstep:-hmax); collect(-hmax:hstep:hmax)]
    len = length(h)
    println(len)
    m = zeros(len)
    
    for i=1:len
        apply_field!(p, h[i])
        m[i] = sum(cos.(p.ϕs))/p.n
    end
    plot(h, m, xlim = (-p.hsat, p.hsat), ylim = (-2,2));
end

function get_magnetization(p::UsovParticle)
    return sum(cos.(p.ϕs))/p.n
end

function magnetize_particle(p::UsovParticle, h::Array{Float64, 1})
    len = length(h)
    
    m = zeros(len)
    for i = 1:len
        apply_field!(p, h[i])
        m[i] = get_magnetization(p)
    end
    
    return m
end