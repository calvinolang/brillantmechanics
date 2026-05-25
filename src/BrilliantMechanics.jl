module BrilliantMechanics

export gravitational_force, kinetic_energy

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

end # module BrilliantMechanics
