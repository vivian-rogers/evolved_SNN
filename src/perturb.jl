#using PlotStuff
using LinearAlgebra
using Random
using Distributions
using Plots
pyplot()


⊗(x,y) = kron(x,y)

function mkfolder(path)
	exists = isdir(path)
	if(exists)
		println("$path already exists...")
		#rm(path, recursive=true)
	else
		mkdir(path)
	end
	return exists
end


function randWeights(P,s)
	A = rand(s,s)
	for i = 1:s
		for j = 1:s
			if(A[i,j] < P)
				A[i,j] = 1
			else
				A[i,j] = 0
			end
		end
	end
	return A
end

function randbin(lo,hi)
	if(rand() > 0.5)
		return lo
	else
		return hi
	end
end

function nnHopping(s)
	A = zeros(s,s)
	inds = [i for i = 1:s]
	inds = shuffle(inds)
	for i = 1:s
		nnR = mod(i+1 - 1,s) + 1
		nnL = mod(i -1 - 1,s) + 1
		A[nnR,i] = randbin(0.5,1)
		A[nnL,i] = randbin(0.5,1)
	end
	return A
end

function randPermMat(s)
	A = zeros(s,s)
	inds = [i for i = 1:s]
	inds = shuffle(inds)
	for i = 1:s
		A[i,inds[i]] = 1
	end
	return A
end

function TMR(ΔV,TMR0)
	γ = 0.5
	return TMR0*γ^2/(ΔV^2 + γ^2)
end

function TMR0grid(σ_TMR,μ_TMR,s)
	TMRs = zeros(s,s)
	TMRdist = Normal(μ_TMR,σ_TMR)
	for i=1:s, j=1:s
		TMRs[i,j] = rand(TMRdist)
	end
	return TMRs
end

function Ggrid(TMR0grid, weightgrid, s, ΔV_in, ΔV_out, μ_TMR, G₁)
	G = zeros(s,s)
	for i=1:s, j=1:s
		TMRᵢⱼ = TMR0grid[i,j]
		ΔV = ΔV_in[i] - ΔV_out[j]
		w = weightgrid[i,j] 
		G[i,j] = w*G₁ + (1-w)*(G₁/(1 + TMR(ΔV,TMRᵢⱼ)))
	end
	return (1/√(s))*G
end



function genΔVs()
	n = 80
	xup = collect(range(0,1,n))
	xdown = collect(range(1,0,n))
	#V_l(x) = x -> 0.5*x
	#V_h(x) = x -> 0.7*x
	#V_0(x) = x -> 0*x
	h = 0.9; l = -0.5
	ΔV_in  = append!(h .* xup,h .* xdown, 0 .* xup, 0 .* xdown, h .* xup,h .* xdown)
	ΔV_out = append!(0 .* xup,0 .* xdown, l .* xup, l .* xdown, l .* xup,l .* xdown)
	return [[ΔV_in[i], ΔV_out[i]] for i = eachindex(ΔV_in) ]
end

function g(x,μ,σ)
	return exp( -(1/2)*((x - μ)/σ)^2 )
end

function ΔVgrid(V_in,V_out,s)
	A = zeros(s,s)
	for i = 1:s, j = 1:s
		A[i,j] = V_in[i] - V_out[j]
	end
	return A
end

function main(name)
	TMR = 1.5; n = 1; s = 30; P = 3/s; G₀ = 1.0; G₁ = G₀*(TMR + 1); U = 0
	path = "./eigval_voltage_sweep"
	mkfolder(path)
	println("Calculating eigvls of reservoir $name, TMR = $TMR, s = $s")
	ΔV₀ = genΔVs()
	TMR0s = TMR0grid(0.0,TMR,s)
	w1 = randWeights(P,Int(s/2))
	w2 = randPermMat(Int(s/2))
	weightgrid = [1 0; 0 0]⊗w2 .+ [0 0; 0 1]⊗w1
	#weightgrid = nnHopping(s)
	#show(weightgrid)
	varr1 = rand(s); varr2 = rand(s);
	hubbardU = 2*U*Diagonal(rand(s) .- 0.5) 
	anim = @animate for ΔV ∈ ΔV₀
		#ΔV_in  = ΔV[1].*[mod(i,2) for i = 1:s]
		#ΔV_out = ΔV[2].*[mod(i+1,2) for i = 1:s]
		ΔV_in = ΔV[1].*[g(i,s/2,√(s)) for i = 1:s]
		ΔV_out = ΔV[2].*[g(i,1,√(s)) for i = 1:s]
		#ΔV_in = ΔV[1].*varr1
		#ΔV_out = ΔV[2].*varr2
		G = Ggrid(TMR0s, weightgrid, s, ΔV_in, ΔV_out, TMR, G₁) .+ hubbardU
		E = (√(s)/2.9)*eigvals(G)
		scatter(real.(E),imag.(E), xlabel = "Re(Eᵢ)", ylabel = "Im(Eᵢ)", size=(600,600))
		xlims!((-1.2,1.2))
		ylims!((-1.2,1.2))
		
		x = y = 1:s; 
		#x = [i for i = 1:s]; y = [i for i = 1:s]; 
		ΔV(i,j) = abs(ΔV_in[i] - ΔV_out[j]); 

		#heatmap!(x,y,G, title="|Gᵢⱼ|", xlims=(1,s), ylims=(1,s), clim=(0,2), inset = (1, bbox(0.05,0.05,0.25,0.25, :top, :right)), subplot = 2)
		heatmap!(x,y,ΔV, title="|ΔVᵢⱼ|", xlims=(1,s), ylims=(1,s), clim=(0,1.5), inset = (1, bbox(0.05,0.05,0.25,0.25, :top, :left)), subplot = 2)
	end
	println("Done animating $name, saving")
	#display(anim)
	#gif(anim,name * ".mp4", fps = 25)
	gif(anim,path*"/"*name * ".gif", fps = 30)
end


main("ex1rand")
main("ex2rand")
main("ex3rand")
main("ex4rand")

# G₀ * TMR = G₁ - G₀
# G₀(TMR + 1)  =  G₁
# G₁*(TMR + 1)^-1 = G₀
#=eigspre = []
eigspost = []
for i = 1:n
	A₀ = (1/C)*randArray(P,s,G₀,G₁)
	E = eigvals(A₀)

	append!(eigspre,[i for i in E])
end
#fig = histogram2d(collect(eigs), nbins = (n, n/2), show_empty_bins = true, normed = true, aspect_ratio = 1)
fig = scatter!(real.(eigs),imag.(eigs))
display(fig) =#
