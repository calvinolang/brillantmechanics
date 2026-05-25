module BrilliantMechanics

export gravitational_force, kinetic_energy, time_to_reach_distance, projectile_range, required_launch_speed, minimum_valid_angle

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

Calculate the first positive time ``t > 0`` (or non-negative time if stationary) at which an object under constant acceleration ``a`` 
and initial velocity ``v_0`` reaches a given distance (displacement) ``d``:

``d = v_0 t + \\frac{1}{2} a t^2``

This solves the quadratic kinematics equation by finding the smallest positive real root. 
If the target is the starting point (``d = 0``) and the object is moving or accelerating, 
it calculates the time the object takes to return to the starting point.

# Arguments
* `d::Real`: The target distance (displacement) in meters (m).
* `v0::Real`: The initial velocity in meters per second (m/s).
* `a::Real`: The constant acceleration in meters per second squared (m/s²).

# Returns
* The first positive time ``t > 0`` (or ``0.0`` if remaining stationary) in seconds (s).

# Errors
* `DomainError`: If the target distance is unreachable, or if both roots are negative (meaning the target was only reached in the past).
"""
function time_to_reach_distance(d::Real, v0::Real, a::Real)
    # If already at the target and not moving/accelerating, time is 0
    if d == 0 && v0 == 0 && a == 0
        return 0.0
    end
    
    if a == 0
        if v0 == 0
            throw(DomainError((d, v0, a), "Target distance is unreachable because velocity and acceleration are both zero."))
        end
        t = d / v0
        # If starting at 0 and moving at constant velocity, we never return to 0 in positive time.
        if t <= 0
            throw(DomainError(t, "Target distance is unreachable in positive future time."))
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
    
    # Filter for strictly positive roots to find future arrival/return times.
    # We use a tiny tolerance (1e-12) to avoid returning t = 0 (the starting point).
    roots = Float64[]
    if t1 > 1e-12
        push!(roots, t1)
    end
    if t2 > 1e-12
        push!(roots, t2)
    end
    
    if isempty(roots)
        # If we are already at the target distance, and there is no future return path, return 0.0
        if d == 0
            return 0.0
        end
        throw(DomainError((t1, t2), "Target distance is unreachable in positive future time."))
    end
    
    return minimum(roots)
end

"""
    projectile_range(a::Real, v0::Real, theta::Real)

Calculate the horizontal distance (range) traveled by a projectile launched from ground level 
before hitting the ground, assuming no air resistance:

``R = \\frac{v_0^2 \\sin(2\\theta)}{g}``

Where ``g = |a|`` is the magnitude of the gravitational acceleration.

# Arguments
* `a::Real`: The gravitational acceleration in meters per second squared (m/s²). The magnitude ``|a|`` is used.
* `v0::Real`: The initial launch speed in meters per second (m/s).
* `theta::Real`: The launch angle relative to the horizontal, in radians (rad).

# Returns
* The horizontal distance (range) traveled in meters (m). If the projectile is launched downward (``\\theta \\leq 0``) or speed is zero, the range is ``0.0``.

# Errors
* `DomainError`: If the acceleration ``a`` is zero (since gravity is required for the projectile to return to the ground).
"""
function projectile_range(a::Real, v0::Real, theta::Real)
    g = abs(a)
    if g == 0
        throw(DomainError(a, "Gravitational acceleration magnitude must be greater than zero."))
    end
    if v0 < 0
        throw(DomainError(v0, "Initial speed v0 must be non-negative."))
    end
    
    # If launch angle points downwards or initial speed is 0, it hits the ground immediately.
    if sin(theta) <= 0 || v0 == 0
        return 0.0
    end
    
    # Range formula: R = (v0^2 * sin(2*theta)) / g
    return Float64((v0^2 * sin(2 * theta)) / g)
end

"""
    required_launch_speed(a::Real, y0::Real, yf::Real, xf::Real, theta::Real)

Calculate the initial launch speed ``v_0`` required for a projectile to reach a target 
final position ``(x_f, y_f)`` when launched from ``(0, y_0)`` at an angle ``\\theta`` 
relative to the horizontal, under constant gravitational acceleration:

``v_0 = \\frac{x_f}{\\cos(\\theta)} \\sqrt{\\frac{g}{2 (y_0 - y_f + x_f \\tan(\\theta))}}``

Where ``g = |a|`` is the magnitude of the gravitational acceleration.

# Arguments
* `a::Real`: The gravitational acceleration in meters per second squared (m/s²). The magnitude ``|a|`` is used.
* `y0::Real`: The initial vertical position (height) in meters (m).
* `yf::Real`: The final target vertical position (height) in meters (m).
* `xf::Real`: The target horizontal distance in meters (m). Must be non-zero.
* `theta::Real`: The launch angle relative to the horizontal, in radians (rad).

# Returns
* The required initial launch speed ``v_0`` in meters per second (m/s).

# Errors
* `DomainError`: If gravity ``a`` is zero.
* `DomainError`: If target horizontal distance ``x_f`` is zero.
* `DomainError`: If the launch angle ``\\theta`` is vertical (i.e. ``\\cos(\\theta) = 0``), which makes horizontal travel impossible.
* `DomainError`: If the straight-line launch path does not pass above the target height (i.e. ``y_0 - y_f + x_f \\tan(\\theta) \\leq 0``), which makes the target physically unreachable under downward gravity.
"""
function required_launch_speed(a::Real, y0::Real, yf::Real, xf::Real, theta::Real)
    g = abs(a)
    if g == 0
        throw(DomainError(a, "Gravitational acceleration magnitude must be greater than zero."))
    end
    if xf == 0
        throw(DomainError(xf, "Target horizontal distance xf must be non-zero."))
    end
    
    cos_theta = cos(theta)
    if abs(cos_theta) < 1e-12
        throw(DomainError(theta, "Launch angle theta cannot be vertical (cos(theta) = 0) when horizontal travel is required."))
    end
    
    # Term in the denominator of the square root
    denom = y0 - yf + xf * tan(theta)
    if denom <= 0
        throw(DomainError(denom, "Target is physically unreachable: the straight-line launch trajectory must pass above the target height yf under downward gravity."))
    end
    
    val = g / (2 * denom)
    return Float64((abs(xf) / abs(cos_theta)) * sqrt(val))
end

"""
    minimum_valid_angle(a::Real, v0::Real, y0::Real, yf::Real, xf::Real, theta_min::Real)

Calculate the minimum launch angle ``\\theta \\geq \\theta_{\\min}`` required for a projectile 
with initial speed ``v_0`` launched from ``(0, y_0)`` to hit a target landing coordinate ``(x_f, y_f)`` 
under gravitational acceleration, satisfying the minimum angle constraint:

This is solved by finding the roots of the quadratic trajectory equation:

``\\frac{g x_f^2}{2 v_0^2} \\tan^2(\\theta) - x_f \\tan(\\theta) + \\left( y_f - y_0 + \\frac{g x_f^2}{2 v_0^2} \\right) = 0``

And then selecting the smallest valid angle ``\\theta \\geq \\theta_{\\min}`` of the two trajectory solutions (high-ball and low-ball).

# Arguments
* `a::Real`: The gravitational acceleration in meters per second squared (m/s²). The magnitude ``|a|`` is used.
* `v0::Real`: The initial launch speed in meters per second (m/s). Must be positive.
* `y0::Real`: The initial vertical position (height) in meters (m).
* `yf::Real`: The target landing vertical position (height) in meters (m).
* `xf::Real`: The target horizontal distance in meters (m). Must be non-zero.
* `theta_min::Real`: The minimum launch angle constraint in radians (rad).

# Returns
* The minimum valid launch angle ``\\theta \\geq \\theta_{\\min}`` in radians (rad).

# Errors
* `DomainError`: If gravity ``a`` is zero.
* `DomainError`: If speed ``v_0`` is non-positive.
* `DomainError`: If target horizontal distance ``x_f`` is zero.
* `DomainError`: If the target coordinate ``(x_f, y_f)`` is physically unreachable with speed ``v_0`` (i.e. outside the parabola of safety).
* `DomainError`: If the minimum angle constraint ``\\theta_{\\min}`` is greater than both physically valid launch angles.
"""
function minimum_valid_angle(a::Real, v0::Real, y0::Real, yf::Real, xf::Real, theta_min::Real)
    g = abs(a)
    if g == 0
        throw(DomainError(a, "Gravitational acceleration magnitude must be greater than zero."))
    end
    if v0 <= 0
        throw(DomainError(v0, "Initial speed v0 must be positive."))
    end
    if xf == 0
        throw(DomainError(xf, "Target horizontal distance xf must be non-zero."))
    end
    
    # Quadratic coefficients for tan(theta): A * z^2 + B * z + C = 0
    A = (g * xf^2) / (2 * v0^2)
    B = -xf
    C = yf - y0 + A
    
    discriminant = B^2 - 4 * A * C
    if discriminant < 0
        throw(DomainError(discriminant, "Target is physically unreachable: the initial speed v0 is too small to reach (xf, yf)."))
    end
    
    sqrt_D = sqrt(discriminant)
    z1 = (-B + sqrt_D) / (2 * A)
    z2 = (-B - sqrt_D) / (2 * A)
    
    # Convert tan(theta) roots to angles. 
    # Since the projectile must travel toward xf, the direction of cos(theta) must match the sign of xf.
    if xf > 0
        theta_a = atan(z1)
        theta_b = atan(z2)
    else
        theta_a = atan(z1) + pi
        theta_b = atan(z2) + pi
    end
    
    # Sort the two launch angles (low-ball and high-ball)
    theta_low = min(theta_a, theta_b)
    theta_high = max(theta_a, theta_b)
    
    if theta_min <= theta_low
        return Float64(theta_low)
    elseif theta_min <= theta_high
        return Float64(theta_high)
    else
        throw(DomainError(theta_min, "The minimum angle constraint theta_min is too large: no valid trajectory angle >= theta_min can hit the target (xf, yf)."))
    end
end

end # module BrilliantMechanics
