using Pkg
Pkg.add("Enzyme")
Pkg.add("Oceananigans")
Pkg.add("UnicodePlots")
Pkg.add("KernelAbstrations")

using Enzyme
using Oceananigans
using UnicodePlots

using KernelAbstractions: @index, @kernel, Event, MultiEvent, NoneEvent
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Utils: launch!
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
        # fld[-Hx+1:0    ,1,1] .= fld[Nx-Hx+1:Nx    ,1,1]
        # fld[ Nx+1:Nx+Hx,1,1] .= fld[      1:1+Hx-1,1,1]
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
  ### for i in 1:nsteps
  ###   del21d!(d2buf,fld)
  ###   fld[1:fld.grid.Nx,1,1].=fld[1:fld.grid.Nx,1,1].+k.*d2buf.*dt
  ###   halo1d!(fld)
  ### end
  ## halo1d!(c)
  for i in 1:nsteps
    ### syntax below fails
    ### fld[10:11,1,1] .= fld[10,1,1]+fld[11,1,1]
    del21d!(d2buf,fld)
    ### kernel style
    event=launch!(arch,grid,:xyz,del21d_k!,d2buf_k,fld)
    wait(device(arch), event) 
    for j in 1:fld.grid.Nx
     d2buf[i] = d2buf_k[i]
    end
    for j in 1:fld.grid.Nx
      fld[j,1,1] = fld[j,1,1] + k*d2buf[j]*dt
    end
    # tmp = fld[10,1,1]+fld[11,1,1]
    # fld[10,1,1] = tmp
    # fld[11,1,1] = tmp
    halo1d!(fld)
    # fld[10,1,1] = fld[10,1,1]+fld[11,1,1]
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

g=grid
lineplot(c[1:g.Nx,1,1],ylim=(0,1))
lineplot(c2[1:g.Nx,1,1]-c[1:g.Nx,1,1])

#
# Lets try some enzymeing
#
bc=CenterField(grid)
set!(c,f)
set!(bc,0)
j=[0.]
bj=[1.]
autodiff(diffuse1d_model!,Duplicated(j,bj), Duplicated(c,bc) )

