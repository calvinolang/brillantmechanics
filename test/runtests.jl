using BrilliantMechanics
using Test

@testset "BrilliantMechanics.jl Tests" begin

    @testset "Gravitational Force" begin
        # Standard calculation
        # F = G * (m1 * m2) / r^2
        # m1 = 1e24, m2 = 1e24, r = 1e8
        # F = 6.67430e-11 * (1e48) / 1e16 = 6.67430e21
        @test gravitational_force(1e24, 1e24, 1e8) ≈ 6.67430e21
        
        # Test unit mass and unit distance
        @test gravitational_force(1, 1, 1) == 6.67430e-11
        
        # Exception handling: zero or negative distance
        @test_throws DomainError gravitational_force(1, 1, 0)
        @test_throws DomainError gravitational_force(1, 1, -5)
        
        # Exception handling: negative mass
        @test_throws DomainError gravitational_force(-10, 1, 5)
        @test_throws DomainError gravitational_force(1, -10, 5)
    end

    @testset "Kinetic Energy" begin
        @test kinetic_energy(2.0, 3.0) ≈ 9.0
        @test kinetic_energy(10.0, 0.0) ≈ 0.0
        
        # Exception handling: negative mass
        @test_throws DomainError kinetic_energy(-5, 2)
    end

    @testset "Time to Reach Distance" begin
        # Stationary start case
        @test time_to_reach_distance(0.0, 0.0, 0.0) ≈ 0.0
        
        # Arrow shot up, coming back to earth (returning to d=0 under gravity)
        # v0 = 10.0 m/s, a = -9.8 m/s^2, should return 2 * v0 / |a|
        @test time_to_reach_distance(0.0, 10.0, -9.8) ≈ (2 * 10.0 / 9.8)
        
        # Constant velocity
        @test time_to_reach_distance(100.0, 10.0, 0.0) ≈ 10.0
        
        # Constant acceleration from rest
        @test time_to_reach_distance(10.0, 0.0, 2.0) ≈ sqrt(10.0)
        
        # Deceleration (two positive roots, should return the first one)
        @test time_to_reach_distance(2.0, 3.0, -2.0) ≈ 1.0
        
        # Unreachable cases (discriminant < 0)
        @test_throws DomainError time_to_reach_distance(10.0, 2.0, -2.0)
        
        # Unreachable cases (zero velocity and acceleration, non-zero distance)
        @test_throws DomainError time_to_reach_distance(10.0, 0.0, 0.0)
        
        # Unreachable cases (negative time required)
        @test_throws DomainError time_to_reach_distance(10.0, -5.0, 0.0)
        # Fallback return case (object starts at 0 and accelerates away, never returning in future, so fallback is 0.0)
        @test time_to_reach_distance(0.0, 5.0, 2.0) ≈ 0.0
    end

    @testset "Projectile Range" begin
        # Standard launch at 45 degrees (pi/4) - yields maximum range
        @test projectile_range(10.0, 10.0, pi/4) ≈ 10.0
        
        # Test unit-agnostic / negative gravity magnitude handling
        @test projectile_range(-9.8, 10.0, pi/4) ≈ (100.0 / 9.8)
        
        # Angle of 0 (fired horizontally along the ground)
        @test projectile_range(9.8, 15.0, 0.0) ≈ 0.0
        
        # Fired downward (theta < 0)
        @test projectile_range(9.8, 15.0, -0.5) ≈ 0.0
        
        # Fired straight up (theta = pi/2) - range should be 0 (re-lands at same spot)
        @test projectile_range(9.8, 15.0, pi/2) ≈ 0.0 atol=1e-12
        
        # Zero speed
        @test projectile_range(9.8, 0.0, pi/4) ≈ 0.0
        
        # Domain errors
        @test_throws DomainError projectile_range(0.0, 10.0, pi/4)  # Zero gravity
        @test_throws DomainError projectile_range(9.8, -5.0, pi/4)  # Negative speed
    end

    @testset "Required Launch Speed" begin
        # Flat ground: v0 = sqrt(g * xf / sin(2*theta))
        # g = 9.81, xf = 10.0, theta = pi/4 (45 deg) -> v0 = sqrt(98.1)
        @test required_launch_speed(9.81, 0.0, 0.0, 10.0, pi/4) ≈ sqrt(98.1)
        
        # Horizontal firing from cliff (theta = 0): v0 = xf * sqrt(g / (2 * y0))
        # g = 10.0, y0 = 10.0, yf = 0.0, xf = 20.0, theta = 0.0 -> v0 = 20.0 * sqrt(10.0 / 20.0) = 10 * sqrt(2)
        @test required_launch_speed(10.0, 10.0, 0.0, 20.0, 0.0) ≈ 10.0 * sqrt(2.0)
        
        # Negative gravity magnitude handling
        @test required_launch_speed(-10.0, 10.0, 0.0, 20.0, 0.0) ≈ 10.0 * sqrt(2.0)
        
        # Fired backward (xf = -20, theta = pi) -> cos(pi) = -1, tan(pi) = 0
        @test required_launch_speed(10.0, 10.0, 0.0, -20.0, pi) ≈ 10.0 * sqrt(2.0)
        
        # Domain error: vertical launch angle (cos(theta) = 0)
        @test_throws DomainError required_launch_speed(9.8, 0.0, 10.0, 10.0, pi/2)
        
        # Domain error: zero gravity
        @test_throws DomainError required_launch_speed(0.0, 0.0, 0.0, 10.0, pi/4)
        
        # Domain error: zero horizontal distance
        @test_throws DomainError required_launch_speed(9.8, 0.0, 0.0, 0.0, pi/4)
        
        # Domain error: physically unreachable (straight line points below landing height)
        # y0 = 0.0, yf = 10.0, xf = 20.0, theta = -pi/4 (tan(theta) = -1) -> denom = 0 - 10 - 20 = -30 <= 0
        @test_throws DomainError required_launch_speed(9.8, 0.0, 10.0, 20.0, -pi/4)
    end

end
