### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ c27858ac-c842-11ea-367b-e17f637879d5
using VoronoiFVM, ExtendableGrids, PyPlot

# ╔═╡ 32ec0b8e-f83a-11ea-15ab-eda28374b333
@bind add_packages html"""Add packages ? <input type=button>"""

# ╔═╡ 24eaafb6-f39b-11ea-20e7-b3728d87f54a
begin
	add_packages 
	using Pkg ; 
	Pkg.add("VoronoiFVM")
	Pkg.add("ExtendableGrids")
	Pkg.add("PlutoUI")
	Pkg.add("PyPlot")
end

# ╔═╡ 2fc26784-c840-11ea-0e74-9dc1d67db495
md"""
### Handling of hetero boundaries

We solve he problem
$ -\Delta d(\mu,E(x))=0$ in $\Omega=(0,1)$ with Dirichlet boundary conditions.
"""

# ╔═╡ d0171e64-f842-11ea-2fee-735051b201e5
md"""
Define the density function $d(\mu,E)$
"""

# ╔═╡ b061020c-f83c-11ea-0d50-9770bb9dfd5e
d(μ,E)=exp(μ/E);

# ╔═╡ e1dfb200-f842-11ea-0343-ffc7780dfe8e
md"""
Set the "reference energy" $E(x)$ such that E(x)=E[1] for x<0.5 and E(x)=E[2] for x>0.5
"""

# ╔═╡ fa643f4a-f83c-11ea-2702-45e5a822d9fe
E=[1.0,1]

# ╔═╡ aabfbe02-f84a-11ea-19a4-2d97430dec2c
md"""
Set the boundary values
"""

# ╔═╡ bbabd2b0-f844-11ea-021b-776ec927437b
bc=[1,2]

# ╔═╡ dbf49d60-f847-11ea-3924-df3513a9799a
html"""<hr>"""

# ╔═╡ 2ecc0eb6-f842-11ea-19fb-7f20249757f8
 md"""
#### Run the solver

Set the number of grid points (it must be even)
"""

# ╔═╡ e793e708-f841-11ea-169f-c908c4eb30a3
n=50

# ╔═╡ f654d6c2-c843-11ea-2889-ab89f8c07f0b
begin
	h=1.0/n
X=collect(0:h:1); 
global grid=simplexgrid(X)
ExtendableGrids.bfacemask!(grid,[0.5],[0.5],3)
ExtendableGrids.cellmask!(grid,[0.0],[0.5],1)
ExtendableGrids.cellmask!(grid,[0.5],[1.0],2)
	grid
end

# ╔═╡ 87267fb2-f3ba-11ea-2422-af9bf7ca5e14
begin
	fig1=PyPlot.figure(figsize=(3,3))
	ExtendableGrids.plot(grid,Plotter=PyPlot)
	fig1
end

# ╔═╡ d82b4036-c844-11ea-14c2-a57826337a32
function g!(f,u,edge)
        f[1]=d(u[1,1],E[edge.region]) - d(u[1,2],E[edge.region])
end;

# ╔═╡ fca38020-c844-11ea-0f9f-e9df10f5cd56
physics=VoronoiFVM.Physics(num_species=1,flux=g!)

# ╔═╡ 1bb53c94-c845-11ea-22f5-dd695a94c4f4
sys=VoronoiFVM.DenseSystem(grid,physics)

# ╔═╡ 16bbef44-c845-11ea-1849-4b9553f89757
enable_species!(sys,1,[1,2])

# ╔═╡ 015b6192-c846-11ea-07c0-f378598ec600
boundary_dirichlet!(sys,1,1,bc[1])

# ╔═╡ 1fd7001c-f83d-11ea-05e7-27edf61f3a93
boundary_dirichlet!(sys,1,2,bc[2])

# ╔═╡ 264d2bb6-c846-11ea-12df-afc41a92fd8b
inival=unknowns(sys); inival.=0.5*(bc[1]+bc[2]);

# ╔═╡ 6cd66908-c846-11ea-0893-658e55bd7058
E;bc;solution0=unknowns(sys);

# ╔═╡ 597b1398-f848-11ea-0783-09022cb7eb1e
solution=solve!(solution0,inival,sys)

# ╔═╡ 870a884e-f849-11ea-055c-d390fcb7a951
begin
	# trigger reactivity
	fig=PyPlot.figure()
	PyPlot.plot(X,solution[1,:],label="μ")
	PyPlot.title("h=$(h), E=$(E)")
	PyPlot.xlabel("x")
	cellnodes=grid[CellNodes]
	cellregions=grid[CellRegions]
	local dmax=0.0
	for icell=1:size(cellnodes,2)
		i1=cellnodes[1,icell]
		i2=cellnodes[2,icell]
		r=cellregions[icell]
		d1=d(solution[1,i1],E[r])
		d2=d(solution[1,i2],E[r])
		label=""
		if icell==1
			label="d"
		end
		PyPlot.plot([X[i1],X[i2]],[d1,d2],"r",label=label)
		dmax=max(dmax,d1,d2)
	end
	PyPlot.plot([0.5,0.5],[0.0,dmax],"k--")
    PyPlot.legend(loc="center left")
	PyPlot.grid()
	
	fig
end

# ╔═╡ Cell order:
# ╟─32ec0b8e-f83a-11ea-15ab-eda28374b333
# ╟─24eaafb6-f39b-11ea-20e7-b3728d87f54a
# ╠═c27858ac-c842-11ea-367b-e17f637879d5
# ╟─2fc26784-c840-11ea-0e74-9dc1d67db495
# ╟─870a884e-f849-11ea-055c-d390fcb7a951
# ╟─d0171e64-f842-11ea-2fee-735051b201e5
# ╠═b061020c-f83c-11ea-0d50-9770bb9dfd5e
# ╟─e1dfb200-f842-11ea-0343-ffc7780dfe8e
# ╠═fa643f4a-f83c-11ea-2702-45e5a822d9fe
# ╟─aabfbe02-f84a-11ea-19a4-2d97430dec2c
# ╠═bbabd2b0-f844-11ea-021b-776ec927437b
# ╟─dbf49d60-f847-11ea-3924-df3513a9799a
# ╟─2ecc0eb6-f842-11ea-19fb-7f20249757f8
# ╠═e793e708-f841-11ea-169f-c908c4eb30a3
# ╠═f654d6c2-c843-11ea-2889-ab89f8c07f0b
# ╟─87267fb2-f3ba-11ea-2422-af9bf7ca5e14
# ╠═d82b4036-c844-11ea-14c2-a57826337a32
# ╠═fca38020-c844-11ea-0f9f-e9df10f5cd56
# ╠═1bb53c94-c845-11ea-22f5-dd695a94c4f4
# ╠═16bbef44-c845-11ea-1849-4b9553f89757
# ╠═015b6192-c846-11ea-07c0-f378598ec600
# ╠═1fd7001c-f83d-11ea-05e7-27edf61f3a93
# ╠═264d2bb6-c846-11ea-12df-afc41a92fd8b
# ╠═6cd66908-c846-11ea-0893-658e55bd7058
# ╠═597b1398-f848-11ea-0783-09022cb7eb1e
