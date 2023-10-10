using JuMP, ExcelReaders, CSV, DataFrames, Distributions, Gurobi, Ipopt, RCall

import XLSX

function evpi(N,scenarios,data1,data2,data3,line,demand,Safl)

    c = 0.75
    p_min = 2000
    p_max = 10000

    ben = 881

    alfa = zeros((N,N))
    alfa[1,:] = (ben/11)*[0 11 17 20]
    alfa[2,:] = (ben/11)*[11 0 22 33]
    alfa[3,:] = (ben/11)*[17 22 0 18]
    alfa[4,:] = (ben/11)*[20 33 18 0]

    Xmax = 1000
    Dmax = demand*Xmax
    Y = 0.25*Safl*Xmax

    x0 = data2[line,:]*Xmax
    d0 = data1[line,:]*Dmax

    valores = zeros((1,100))
    opt_gap = zeros((1,100))
    precios = zeros((N,100))
    concentraciones = zeros((N,100))
    flujos = zeros((N^2,100))
    dem_eff = zeros((N,100))
    infactible = zeros((1,100))

    t3 = time_ns()/1.0e9

    for sim=1:10
        y = data3[sim,1:4]*Y
        muest = data3[sim,5:8]

        #matriz de escenarios, por generalizar
        D0 = zeros((N,scenarios))

        for i=1:N
            for k=1:scenarios
                D0[i,k] = (0.7+0.3*(k-1))*d0[i]
            end
        end

        delta = 0.9

        m = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "OutputFlag"=>0, "TimeLimit" => 1800))

        M = 4*p_max
        K = p_max*(Xmax + 1.3*Dmax + Y)

        @variable(m, x[1:N] >=0)
        @variable(m, v[1:N,1:N] >=0, Int)
        @variable(m, p_min <= p[1:N] <= p_max)

        @variable(m, gamma[1:N])
        @variable(m, lambda[1:N,1:N])

        @variable(m, -1 <= beta[1:N,1:scenarios] <= 1)

        #desacople HC gamma y lambda
        @variable(m, 0 <= w[1:N] <= 1, Int)
        @variable(m, 0 <= z[1:N,1:N] <= 1, Int)

        #linealizacion por pedazos función objetivo
        @variable(m, 0 <= a[1:N] <= 1, Int)
        @variable(m, 0 <= t1[1:N] <= K)
        @variable(m, 0 <= t2[1:N] <= K)

        #linealización por pedazos betas
        @variable(m, 0 <= b1[1:N,1:scenarios] <= 1, Int)
        @variable(m, 0 <= b2[1:N,1:scenarios] <= 1, Int)
        @variable(m, 0 <= b3[1:N,1:scenarios] <= 1, Int)
        @variable(m, -K <= s1[1:N,1:scenarios] <= K)
        @variable(m, -K <= s2[1:N,1:scenarios] <= K)
        @variable(m, -K <= s3[1:N,1:scenarios] <= K)

        @expression(m, dEff[i=1:N], muest[i]*d0[i]*(1-delta*((p[i]-p_min)/(p_max-p_min))))
        @expression(m, dEff_matrix[i=1:N,k=1:scenarios], D0[i,k]*(1-delta*((p[i]-p_min)/(p_max-p_min))))

        #Zonas: Santiago (1), Renca (2), Maipú (3), La Florida (4), Las Condes (5)
        #\alpha se calcula con distancia entre "centros" y tomando performance promedio: 11Km/Litro
        #Precio Promedio Bencina 95 en RM: $881 (20-05-2021)
        #Las distancias no son necesariamente simétricas, ¿considerar?

        @objective(m, Max, (1-c)*sum(p[i]*(-K*a[i] + t1[i]) + p[i]*dEff[i] for i=1:N))

        @constraint(m, [i=1:N],  t1[i] - K*a[i] <= 0)
        @constraint(m, [i=1:N],  t2[i] - K*(1 - a[i]) <= 0)

        @constraint(m, [i=1:N],  -K*a[i] + t1[i] + t2[i] == x[i] + y[i] - dEff[i])

        #Definición variable x
        @constraint(m, [i=1:N],  x[i] == x0[i] + sum(v[j,i] - v[i,j] for j=1:N))
        @constraint(m, [i=1:N],  v[i,i] == 0)

        #######################################################

        #Ecuaciones KKT_nodo1

        for i=1:N
            for j=1:N
                if(i==j)
                    continue
                end
                #KKTs por nodo
                @constraint(m, sum(p[i]*beta[i,k]-p[j]*beta[j,k] for k=1:scenarios) + alfa[i,j] - lambda[i,j] + gamma[i] == 0)
                #Desacople HC lambda
                @constraint(m, 0 <= M*z[i,j] + lambda[i,j])
                @constraint(m, 0 >= -M*z[i,j] + lambda[i,j])
                #Desacople HC v
                @constraint(m, 0 <= M*(1-z[i,j]) + v[i,j])
                @constraint(m, 0 >= -M*(1-z[i,j]) + v[i,j])

            end
            #Desacople HC gammas
            @constraint(m, 0 <= gamma[i] + M*w[i])
            @constraint(m, 0 >= gamma[i] - M*w[i])
            #Desacople v-x0
            @constraint(m, 0 <= sum(v[i,j] for j=1:N) - x0[i] + M*(1-w[i]))
            @constraint(m, sum(v[i,j] for j=1:N) - x0[i] - M*(1-w[i]) <= 0)

            for k=1:scenarios
                @constraint(m, b1[i,k]+b2[i,k]+b3[i,k] == 1)
                @constraint(m, s1[i,k]+s2[i,k]+s3[i,k] == dEff_matrix[i,k]-x[i])
                @constraint(m, s1[i,k] + K*b1[i,k] >= 0)
                @constraint(m, s1[i,k] <= 0)
                @constraint(m, s2[i,k] >= 0)
                @constraint(m, s2[i,k] - Y*b2[i,k] <= 0)
                @constraint(m, s3[i,k] - Y*b3[i,k] >= 0)
                @constraint(m, s3[i,k] - K*b3[i,k] <= 0)
                @constraint(m, beta[i,k] == -c*(1/scenarios)*(s2[i,k]/Y + b3[i,k]))
            end

        end

        #######################################################

        status = optimize!(m)

        avance = string("EVPI", " Linea ","$line",", Dmax ", "$demand",", Y ","$Safl")

        println(avance)
        println(sim)

        if length(string(JuMP.termination_status(m))) == 10
            infactible[sim] = sim
            continue

        elseif length(string(JuMP.termination_status(m))) != 10
            valores[sim] = JuMP.objective_value(m)
            opt_gap[sim] = JuMP.relative_gap(m)
            precios[1:N,sim] = JuMP.value.(p)
            concentraciones[1:N,sim] = JuMP.value.(x)
            for i=1:N
                flujos[(i-1)*N+1:i*N,sim] = JuMP.value.(v[i,:])
            end

            dem_eff[1:N,sim] = JuMP.value.(dEff)
        end

    end

    t2 = time_ns()/1.0e9
    tf = t2-t3

    datos = zeros((1,2))

    datos[1] = tf
    datos[2] = sum(opt_gap)/10

    Data1 = DataFrame(valores,:auto);
    Data2 = DataFrame(precios,:auto);
    Data3 = DataFrame(concentraciones,:auto);
    Data4 = DataFrame(flujos,:auto);
    Data5 = DataFrame(dem_eff,:auto);
    Data6 = DataFrame(infactible,:auto);

    sufijo = string("line","$line","Dmax","$demand","_","Y","$Safl")

    cd("EV_Perfect")

    CSV.write(string("Resultados_evpi1_",sufijo,".csv"), Data1)
    CSV.write(string("Resultados_evpi2_",sufijo,".csv"), Data2)
    CSV.write(string("Resultados_evpi3_",sufijo,".csv"), Data3)
    CSV.write(string("Resultados_evpi4_",sufijo,".csv"), Data4)
    CSV.write(string("Resultados_evpi5_",sufijo,".csv"), Data5)
    CSV.write(string("Resultados_evpi6_",sufijo,".csv"), Data6)

    cd("..")

end
