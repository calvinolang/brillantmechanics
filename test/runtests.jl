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
        # Immediate case
        @test time_to_reach_distance(0.0, 5.0, 2.0) ≈ 0.0
        
        # Constant velocity
        @test time_to_reach_distance(100.0, 10.0, 0.0) ≈ 10.0
        
        # Constant acceleration from rest
        @test time_to_reach_distance(10.0, 0.0, 2.0) ≈ sqrt(10.0)
        
        # Deceleration (two positive roots, should return the first one)
        @test time_to_reach_distance(2.0, 3.0, -2.0) ≈ 1.0
        
        # Unreachable cases (discriminant < 0)
        @test_throws DomainError time_to_reach_distance(10.0, 2.0, -2.0)
        
        # Unreachable cases (zero velocity and acceleration)
        @test_throws DomainError time_to_reach_distance(10.0, 0.0, 0.0)
        
        # Unreachable cases (negative time required)
        @test_throws DomainError time_to_reach_distance(10.0, -5.0, 0.0)
    end

end
