using Pkg
Pkg.add(name="Oceananigans",version="0.79.2")
Pkg.add("Plots")
Pkg.add(name="Enzyme",version="0.10.18")
Pkg.add("LaTeXStrings")

using Enzyme
using Oceananigans
using LaTeXStrings
using Plots
# ENV["GKSwstype"]="nul"

using Oceananigans.Architectures: device

# Lets create a one-dimensional array, with N elements, using Oceananigans data structures
#

# Specify computation parameters
arch=CPU()
FT=Float64

# We need a grid to create an Oceananigans array
N = 100
topo = (Periodic, Flat, Flat)
grid = RectilinearGrid(arch, FT, topology=topo, size=(N), halo=2, x=(-1, 1), y=(-1, 1), z=(-1, 1))

# Lets create a couple of non-KA operators and test those first
# 1. del^2
function del21d!(d2buf,fld::Field)
 g  = fld.grid
 Nx = g.Nx
 for i=1:Nx
  d2buf[i] = fld[i-1,1,1]+fld[i+1,1,1]-2fld[i,1,1]
 end
 return
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
        # fld[-Hx+1:0    ,1,1] .= fld[Nx-Hx+1:Nx    ,1,1]
        # fld[ Nx+1:Nx+Hx,1,1] .= fld[      1:1+Hx-1,1,1]
        return nothing
end
# 3. simple model
function diffuse1d_model!(jcost,fld)
  grid=fld.grid
  d2buf=ones(grid.Nx)
  arch=grid.architecture
  k=1.0
  dt=0.25
  nsteps=50
  for i in 1:nsteps
    del21d!(d2buf,fld)
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
# set!(c,0)
# c[15,1,1]=1.
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


g=grid
# lineplot(c[1:g.Nx,1,1],ylim=(0,1))
# lineplot(c2[1:g.Nx,1,1]-c[1:g.Nx,1,1])

#
# Lets try some enzymeing
#
bc=CenterField(grid)
set!(c,f)
halo1d!(c)
set!(bc,0)
halo1d!(bc)
j=[0.]
bj=[1.]
autodiff(diffuse1d_model!,Duplicated(j,bj), Duplicated(c,bc) )

halo1d!(bc)
# Plots.plot(bc.data[:,1,1])
ttl=L"\left[J=\phi(n=N_t,i=i_J), i_J=15, N_t=50\right]"
Plots.plot(bc.data[:,1,1];title=ttl,ylabel=L"\texttt{bc}=\frac{\partial J}{\partial \phi(n=0,i)}",xlabel="i",legend=false)
Plots.savefig("no-ka-pre-changes.png")
Nx=grid.Nx
scatter(bc[1:Nx,1,1];markersize=1,legend=false,title=ttl,ylabel=L"\texttt{bc}=\frac{\partial J}{\partial \phi(n=0,i)}",xlabel="i")
Plots.savefig("no-ka-pre-changes-scatter.png")

# dev form
## autodiff(Reverse, diffuse1d_model!,Duplicated(j,bj), Duplicated(c,bc) )

