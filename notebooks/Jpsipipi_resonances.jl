### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 336daff4-5f29-11ef-3032-733d15dff8a2
# ╠═╡ show_logs = false
begin
	using Pkg
	
	Pkg.add([
		Pkg.PackageSpec("DataFrames"),
		Pkg.PackageSpec("JSON"),
		Pkg.PackageSpec(
		url="https://github.com/mmikhasenko/HadronicLineshapes.jl#v0.4.1")])
	
	using HadronicLineshapes
	using DataFrames
	using JSON
end

# ╔═╡ ba65fe31-1ebf-40eb-bdee-74814df06de7
md"""
## f0(1370)

- parametrization `BW(1.37, 0.35)`
- amplitude (-0.751892,-0.320262)

using same code for any other BW in that fit
```c++
complex<double> amp = complex<double>(1.0,0.0) /
                       complex<double>( mass2 - m_mass*m_mass, m_mass*m_width);
```

## Sigma

- amplitude (-1.24404,0.984936)

using:

```c++
double GammaX= sqrt(1.0-4.0*mpi*mpi/s)*m_width;
complex<double> result = complex<double>(1.0,0.0) /
          (complex<double>(m_mass*m_mass-s, 0.0) - complex<double>(0.0, 1.0)*sqrt(s)*GammaX);
```


## f0(980)

- parameters: `mass 0.965 g1 0.165 g2/g1 4.21`
- amplitude (1.65895,0.463452) 

using
```c++
double arg1 = (s-pow(0.13957061+0.13957061,2))*s;
complex<double> rho_pipi;
if(arg1 > 0)
	rho_pipi = complex<double>(sqrt(arg1)/s,0.);
else
    rho_pipi = complex<double>(0.,sqrt(-arg1)/s);

double arg2 = (s-pow(0.493677+0.493677,2))*s;
complex<double> rho_KK;
if(arg2 > 0)
    rho_KK = complex<double>(sqrt(arg2)/s,0.);
else
    rho_KK = complex<double>(0.,sqrt(-arg2)/s);

complex<double> cgPIPI = complex<double>(m_g1,0.)*rho_pipi;
complex<double> cgKK = complex<double>(m_g1*m_g2_over_g1,0.)*rho_KK;


complex<double> amp = complex<double>(1.,0.) / (complex<double>(pow(m_mass,2)-s,0.) 	- complex<double>(0.,1.)*(cgPIPI + cgKK));
```
"""

# ╔═╡ 5fa81bc5-6964-4b95-99dc-5eb1db11b9f1
md"""
## Data file
"""

# ╔═╡ 526226cc-0243-4dbd-8c49-847a36653510
begin
	refs = [
		"1370" => "(-0.751892,-0.320262)",
		"Sigma" => "(-1.24404,0.984936)",
		"980" => "(1.65895,0.463452)"]
	# 
	const mpipi = 0.86669475698368
	const σ_ref = mpipi^2
end

# ╔═╡ 74cf9d20-33ea-4e9e-80d4-c048383eb5fb
begin
	df = DataFrame(name=first.(refs), vals = last.(refs))
	select!(df, :name, :vals => ByRow() do x
		eval(Meta.parse("complex"*x))
	end => :reference_value)
	df.computed_value .= 0.0im
	df
end

# ╔═╡ b47f00c6-8fe2-4912-a518-8c59e3b6dc4b
md"""
## Implementations
"""

# ╔═╡ 88f669f5-8e26-468a-9e78-667c1e4688c4
X_f01370 = BreitWigner(1.37, 0.35);

# ╔═╡ cd4a0d0a-804b-464e-adad-b6824f051970
X_Sigma(σ) = let m = 0.507, Γ = 0.475, mpi = 0.1349766
	1/(m^2 - σ - 1im*Γ*sqrt(σ)*sqrt(1-4mpi^2/σ))
end

# ╔═╡ 60d13c56-73b6-4733-8b2c-5352af9dadb0
X_BW980 = let
	m = 0.965
	m_g1 = 0.165
	m_g2_over_g1 = 4.21
	mπ = 0.13957061
	mK = 0.493677
	# 
	Flatte(m, 
		m_g1, mπ, mπ,
	    m_g2_over_g1*m_g1, mK, mK)
end

# ╔═╡ 1829c73e-3582-4ec7-807e-660fcb73ce8b
md"""
## Comparison
"""

# ╔═╡ fd6a9934-d0cb-4f11-9953-a010b320c9d7
begin
	df[df.name .== "980", :computed_value] .= X_BW980(σ_ref)
	df[df.name .== "1370", :computed_value] .= X_f01370(σ_ref)
	df[df.name .== "Sigma", :computed_value] .= X_Sigma(σ_ref)
	# 
	transform(df, :name,
			[:reference_value, :computed_value] => ByRow((x,y)->round(y/x; digits=3)) =>:ratio )
end

# ╔═╡ Cell order:
# ╠═336daff4-5f29-11ef-3032-733d15dff8a2
# ╟─ba65fe31-1ebf-40eb-bdee-74814df06de7
# ╟─5fa81bc5-6964-4b95-99dc-5eb1db11b9f1
# ╠═526226cc-0243-4dbd-8c49-847a36653510
# ╠═74cf9d20-33ea-4e9e-80d4-c048383eb5fb
# ╟─b47f00c6-8fe2-4912-a518-8c59e3b6dc4b
# ╠═88f669f5-8e26-468a-9e78-667c1e4688c4
# ╠═cd4a0d0a-804b-464e-adad-b6824f051970
# ╠═60d13c56-73b6-4733-8b2c-5352af9dadb0
# ╟─1829c73e-3582-4ec7-807e-660fcb73ce8b
# ╠═fd6a9934-d0cb-4f11-9953-a010b320c9d7
