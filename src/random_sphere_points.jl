function get_points_spherical_random(n::Int)
    r = ones(n)
    θ = acos.(1-2rand(n))
    ϕ = 2pi*rand(n)
    
    return [r θ ϕ]
end

function get_points_spherical_evenly(n::Int)
    r = ones(n)
    θ = acos.(1-2linspace(0.0, 1.0, n))
    ϕ = 2pi*linspace(0.0, 1.0, n)
    
    return [r θ ϕ]
end

function get_points_cartesian_random(n::Int)
    points = get_points_spherical_random(n)
    
    return [sin.(points[:,2]).*cos.(points[:,3]) sin.(points[:,2]).*sin.(points[:,3]) cos.(points[:,2])]
end

function random_sphere_pack(r::AbstractFloat, R::AbstractFloat, n::Int)
        
    step = r/10
    x = [ collect(-step : -step : -R); collect(0 : step : R)]
    y = [ collect(-step : -step : -R); collect(0 : step : R)]
    z = [ collect(-step : -step : -R); collect(0 : step : R)]
    len = length(x)
    places = zeros(len, len, len)
    
    for i = 1:len
        for j = 1:len
            for k = 1:len
                if sqrt(x[i]^2 + y[j]^2 + z[k]^2) < R
                    places[i,j,k] = 1
                end
            end
        end
    end
    
    points = zeros(3, n)
    count = 0
    
    while count < n
        vacant_points = find(places)
        if isempty(vacant_points)
            error("there are no more vacant points in the space of the cluster at `count`=$count")
        end
        
        random_vacant_point = vacant_points[rand(1:length(vacant_points))]
        (I,J,K) = ind2sub(places, random_vacant_point)
        
        if places[I,J,K] != 1
            display("Something went wrong! places(p(1),p(2),p(3))!=1")
        end
        
        for i = 1:len
            for j = 1:len
                for k = 1:len
                    if sqrt((x[i]-x[I])^2 + (y[j]-y[J])^2 + (z[k]-z[K])^2) < 2r
                        places[i,j,k] = 0
                    end
                end
            end
        end
        
        count += 1;
        points[:, count] = [x[I] y[J] z[K]]
    end
    
    return points
end
