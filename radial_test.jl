#####################################################
#####################################################
using Plots
using Printf
#####################################################
#####################################################
function plotpaths(n,n2)

#get n paths with n2 points plot with domain and rho

#clear current plot so they don't overlay
Plots.CURRENT_PLOT.nullableplot = nothing

f(x,y)=2-(x^2+y^2) #mock radial charge density

p=Matrix{Float64}(undef,n2,2*n) #array to hold gps

htheta=(2*pi)/n #theta spacing between gps

intradius=.05 #radius of circle where gps begin

t=collect(range(intradius,1,n2)) #parameter for gps

# Create a grid of x and y values for plotting rho
xx = range(-1.5, 1.5, length=100)
yy = range(-1.5, 1.5, length=100)

# Evaluate the rho on the grid
z = [f(i, j) for i in xx, j in yy]

# Create the contour plot
fig=contour(xx, yy, z, levels=50, color=:plasma, xlabel="x", ylabel="y", title="domain, f(x,y)=2-x^2+y^2 & gradient paths")

for i=1:n #loop to get gps and plot
x=cos(htheta*(i-1))*t
y=sin(htheta*(i-1))*t
plot!(x,y, linewidth=2,legend=false)
end


tt=collect(range(0,2*pi,50)) #parameter for plotting circles

x=cos.(tt)
y=sin.(tt)
plot!(x,y, linewidth=5,legend=false,c=:black) #plot unit cirlce
x=intradius*cos.(tt)
y=intradius*sin.(tt)
plot!(x,y, linewidth=1,legend=false,c=:black,linestyle=:dash) #plot inner circle

return fig

end

#######################################################
#######################################################
function getpaths(n,n2,intradius)
#get n paths with n2 points radius of circle where paths start is intradius
#assume radial charge density e.g. f(x,y)=2-x^2+y^2

p=Matrix{Float64}(undef,n2,2*n) #array to hold gps

htheta=(2*pi)/n #theta spacing between gps

t=collect(range(intradius,1,n2)) #parameter for gps

for i=1:n #loop to get gps
  x=cos(htheta*(i-1))*t
  y=sin(htheta*(i-1))*t
  p[:,2*i-1]=x
  p[:,2*i]=y
end

return p

end
#######################################################
#######################################################
function int_area(n,n2,intradius)

f(x,y)=2-(x^2+y^2) #mock charge density never actually used in code
curv(x,y)=1/sqrt(x^2+y^2) #curvature of isosurface is reciprocal of distance to origin

exact=pi #exact solution area of circle of radius one
intarea=(pi)*intradius^2 #area of interior circle where gps start

p=getpaths(n,n2,intradius) #get array of gradient paths

arc=Matrix{Float64}(undef,n2,n) #array to hold arc lengths l1(s)

h=(1-intradius)/(n2-1) #spacing in gp points

arc[1,:].=(1/n)*pi*2*intradius #l1(s0)

for i=1:n #loop over gps
for j=2:n2 #loop over points for gp
    int=0
    for k=2:j
    int+=curv(p[k,2*i-1],p[k,2*i])*h #integrate curvature rh rule
    end
    arc[j,i]=arc[1,i]*exp(int) #take exponential and multiply by l1(s0)
end
end

#trapezoidal rule to integrate l1(s) down each path
area=sum(arc[1,:])+sum(arc[end,:]) #first and last terms

for i=2:n2-1
area+=sum(arc[i,:])*2 #middle terms
end

area*=h/2 #weight
area+=intarea #add area of interior circle where gps start
err=abs(area-pi)/pi #calculate relative error
return err
end
#######################################################
#######################################################
#####################################################
function errortable()
@printf "\n Approximate area of unit circle using gp-only integration"
@printf "\n Interior circle radius where gps start: .05"
@printf "\n Number of gradient paths: 10"
@printf "\n n is number of points on each path"
@printf "\n Error is relative error"
@printf "\n EOC is estimated order of convergence"
@printf "\n ----------------------------------------------------------------"
@printf "\n n           rel. error         EOC"

n=(10,20,40,80,160,320,640,1280,2560)
er=Array{Float64}(undef,9)
eoc=Array{Float64}(undef,9)

er[1]=int_area(10,n[1],.05)
@printf "\n %g           %g          n/a" n[1] er[1]

for i=2:9
er[i]=int_area(10,n[i],.05)
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" n[i] er[i] eoc[i]
end

end
