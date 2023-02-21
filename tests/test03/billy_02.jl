using Enzyme
using Oceananigans
using KernelAbstractions:  @index, @kernel
using Oceananigans.Utils: launch!
using Oceananigans.Architectures: device
using Plots
using LaTeXStrings

Enzyme.API.looseTypeAnalysis!(true)
Enzyme.API.printall!(true)
EnzymeRules.inactive(::typeof(Oceananigans.Utils.flatten_reduced_dimensions), x...) = nothing
EnzymeRules.inactive(::typeof(Oceananigans.Grids.default_indices), x...) = nothing

# ENV["GKSwstype"]="nul"

arch=CPU()
FT=Float64

N = 100
topo = (Periodic, Flat, Flat)
grid = RectilinearGrid(arch, FT, topology=topo, size=(N), halo=2, x=(-1, 1), y=(-1, 1), z=(-1, 1))

@kernel function del21d_k!(d2buf_k,fld)
  i,j,k = @index(Global,NTuple)
  @inbounds d2buf_k[i,j,k] = fld[i-1,j,k]+fld[i+1,j,k]-2fld[i,j,k]
end

# 2. halo
function halo1d!(fld::Field)
        g  = fld.grid
        Hx = g.Hx
        Nx = g.Nx
        for i=1:Hx
         fld[i-Hx,1,1]=fld[Nx-Hx+i,1,1]
         fld[Nx+i,1,1]=fld[i      ,1,1]
        end
        return nothing
end
# 3. simple model
function diffuse1d_model!(jcost,fld)
  grid=fld.grid
  d2buf=ones(grid.Nx)
  arch=grid.architecture
  d2buf_k = CenterField(grid)
  k=1.0
  dt=0.1
  nsteps=50
  for i in 1:nsteps
    # i = 1
    # del21d!(d2buf,fld)
    for i=1:fld.grid.Nx
      d2buf[i] = fld[i-1,1,1]+fld[i+1,1,1]-2fld[i,1,1]
    end
    ### kernel style 
    workgroup, worksize = Oceananigans.Utils.work_layout(grid, :xyz;
                                      include_right_boundaries=false,
                                      reduced_dimensions=(),
                                      location=nothing, 
                                      only_active_cells=false)

    loop! = del21d_k!(Oceananigans.Architectures.device(arch), workgroup, worksize)
    event = loop!(d2buf_k, fld) #; dependencies=dependencies)
    
    wait(device(arch), event) 
     
    for j in 1:fld.grid.Nx
      d2buf[i] = d2buf_k[i]
    end
    for j in 1:fld.grid.Nx
      fld[j,1,1] = fld[j,1,1] + k*d2buf[j]*dt
    end
    halo1d!(fld)
  end
  jcost[1]=fld[15,1,1].*fld[15,1,1]
  return nothing
end

c     = CenterField(grid)
c2    = CenterField(grid)
f(x, y, z) = exp( -50((x-grid.xᶠᵃᵃ[1])/grid.Lx-0.5)^2 )
set!(c,f)
halo1d!(c)
set!(c2,f)
halo1d!(c2)
j=[0.]
diffuse1d_model!(j,c)
Nx=grid.Nx
ttl=L"\texttt{c2}=\phi(n=0,i),\quad \texttt{c}=\phi(n=N_t,i)"
scatter(c2[1:Nx,1,1];markersize=2,title=ttl,ylabel=L"\phi",xlabel="i",labels="c2")
scatter!(c[1:Nx,1,1];markersize=2,title=ttl,ylabel=L"\phi",xlabel="i",labels="c")
Plots.savefig("checking-operator-01.png")
ttl=L"\texttt{c2}=\phi(n=0,i),\quad \texttt{c}=\phi(n=N_t,i)"
scatter(c2[1:Nx,1,1].-c[1:Nx,1,1];markersize=2,title=ttl,ylabel=L"\phi^{n=0}-\phi^{n=N_t}",xlabel="i",labels="c2 - c")
Plots.savefig("checking-operator-02.png")


bc=CenterField(grid)
set!(c,f)
set!(bc,0)
j=[0.]
bj=[1.]
autodiff(Reverse, diffuse1d_model!,Duplicated(j,bj), Duplicated(c,bc) )

plot(bc.data[:,1,1])
savefig("using-ka.png")
