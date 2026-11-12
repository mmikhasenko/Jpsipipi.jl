using Test

name = only(ARGS)
include(joinpath(@__DIR__, "..", "notebooks", "Jpsipipi_$(name).jl"))

if name == "amplitudes"
    @testset "Reference helicity amplitudes" begin
        comparisons = vcat(f0_comparison, z_comparison)
        @test size(comparisons, 1) == 11
        @test Set(comparisons.name_ls_LS) == Set(keys(checklist["chains"]["amplitude"]))
        # The supplied angles and complex references are rounded decimal values.
        # Compare complex amplitudes, including phase, at their input precision.
        for row in eachrow(comparisons)
            @test isapprox(row.computed_value, row.reference_value; atol = 5e-6, rtol = 0)
        end
    end
elseif name == "angles"
    @testset "Reference helicity angles" begin
        @test size(angle_comparison, 1) == 9
        @test !any(ismissing, angle_comparison.theta_here)
        @test !any(ismissing, angle_comparison.phi_here)
        # Four-vectors are rounded to six decimal places; theta is in degrees.
        for row in eachrow(angle_comparison)
            @test isapprox(row.theta_here, row.theta_ref; atol = 5e-4, rtol = 0)
            # The f0 topology's psi is the second daughter; see the notebook's
            # frame-convention explanation and decayangle Node.boost particle-2 rule.
            # https://github.com/KaiHabermann/decayangle/blob/v1.1.5/src/decayangle/decay_topology.py
            offset = row.triple_name == "psi --> mup mum : topology_f0" ? π : 0
            @test abs(rem2pi(row.phi_here - row.phi_ref - offset, RoundNearest)) < 5e-5
        end
    end
end
