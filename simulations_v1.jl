using JuMP, ExcelReaders, CSV, DataFrames, Distributions, Gurobi, Ipopt, RCall

import XLSX

c = 0.75
p_min = 2000
p_max = 10000

N = 4
scenarios = 3

f = openxl("simulations_d-x.xls")
g = openxl("samples_xi.xls")

include("evpi.jl")
include("evsi.jl")

data1 = readxlsheet(f,"d0")
data2 = readxlsheet(f,"x0")
data3 = readxlsheet(g,"xi_y_d")

ben = 881

alfa = zeros((N,N))
alfa[1,:] = (ben/11)*[0 11 17 20]
alfa[2,:] = (ben/11)*[11 0 22 33]
alfa[3,:] = (ben/11)*[17 22 0 18]
alfa[4,:] = (ben/11)*[20 33 18 0]


for line=1:100
    for demand=1:5
        for Safl=1:3
            evpi(N,scenarios,data1,data2,data3,line,demand,Safl)
            evsi(N,data1,data2,data3,line,demand,Safl)
        end
    end
end
