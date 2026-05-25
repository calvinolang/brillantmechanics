module BrilliantMechanics

export gravitational_force, kinetic_energy, time_to_reach_distance

# --- Constants ---

"""
    G

The Newtonian constant of gravitation in units of m³/(kg·s²).
Reference value: 6.67430e-11.
"""
const G = 6.67430e-11

# --- Mechanics Functions ---

"""
    gravitational_force(m1, m2, r)

Calculate the magnitude of the gravitational force between two point masses `m1` and `m2` 
separated by a distance `r` using Newton's Law of Universal Gravitation:

``F = G \\frac{m_1 m_2}{r^2}``

# Arguments
* `m1::Real`: Mass of the first object in kilograms (kg).
* `m2::Real`: Mass of the second object in kilograms (kg).
* `r::Real`: Distance between the centers of the masses in meters (m).

# Examples
```julia
# Force between two 1 kg masses separated by 1 meter
gravitational_force(1.0, 1.0, 1.0) # approx 6.6743e-11 N
```
"""
function gravitational_force(m1::Real, m2::Real, r::Real)
    if r <= 0
        throw(DomainError(r, "Distance r must be positive and non-zero."))
    end
    if m1 < 0 || m2 < 0
        throw(DomainError((m1, m2), "Masses must be non-negative."))
    end
    
    return G * (m1 * m2) / r^2
end

"""
    kinetic_energy(m, v)

Calculate the kinetic energy of a point object of mass `m` moving with speed `v`:

``E_k = \\frac{1}{2} m v^2``

# Arguments
* `m::Real`: Mass of the object in kilograms (kg).
* `v::Real`: Speed of the object in meters per second (m/s).

# Examples
```julia
# Kinetic energy of a 2 kg object moving at 3 m/s
kinetic_energy(2.0, 3.0) # 9.0 J
```
"""
function kinetic_energy(m::Real, v::Real)
    if m < 0
        throw(DomainError(m, "Mass must be non-negative."))
    end
    
    return 0.5 * m * v^2
end

"""
    time_to_reach_distance(d::Real, v0::Real, a::Real)

Calculate the first non-negative time ``t`` at which an object under constant acceleration ``a`` 
and initial velocity ``v_0`` reaches a given distance (displacement) ``d``:

``d = v_0 t + \\frac{1}{2} a t^2``

This is solved by finding the smallest non-negative root of the quadratic equation:

``\\frac{1}{2} a t^2 + v_0 t - d = 0``

# Arguments
* `d::Real`: The target distance (displacement) in meters (m).
* `v0::Real`: The initial velocity in meters per second (m/s).
* `a::Real`: The constant acceleration in meters per second squared (m/s²).

# Returns
* The first non-negative time ``t \\geq 0`` in seconds (s) to reach the distance.

# Errors
* `DomainError`: If the target distance is unreachable, or if both roots are negative (meaning the target was only reached in the past).
"""
function time_to_reach_distance(d::Real, v0::Real, a::Real)
    if d == 0
        return 0.0
    end
    
    if a == 0
        if v0 == 0
            throw(DomainError((d, v0, a), "Target distance is unreachable because velocity and acceleration are both zero."))
        end
        t = d / v0
        if t < 0
            throw(DomainError(t, "Target distance is unreachable in future time (requires moving backward)."))
        end
        return Float64(t)
    end
    
    # Quadratic equation: 0.5 * a * t^2 + v0 * t - d = 0
    # Discriminant D = v0^2 - 4 * (0.5 * a) * (-d) = v0^2 + 2 * a * d
    discriminant = v0^2 + 2 * a * d
    
    if discriminant < 0
        throw(DomainError(discriminant, "Target distance is unreachable: the object turns back before reaching d."))
    end
    
    sqrt_D = sqrt(discriminant)
    t1 = (-v0 + sqrt_D) / a
    t2 = (-v0 - sqrt_D) / a
    
    # Filter for non-negative roots
    roots = Float64[]
    if t1 >= 0
        push!(roots, t1)
    end
    if t2 >= 0
        push!(roots, t2)
    end
    
    if isempty(roots)
        throw(DomainError((t1, t2), "Target distance is unreachable in future time (both roots are negative)."))
    end
    
    return minimum(roots)
end

end # module BrilliantMechanics
