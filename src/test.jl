using Plots
gr()
N = 10
x = rand(1:10,N)
y = rand(1:10,N)
u = rand(N)
v = rand(N)
scatter(x,y)
fig = quiver!(x,y,quiver=(u,v))

display(fig)
