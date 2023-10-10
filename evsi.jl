using JuMP, ExcelReaders, CSV, DataFrames, Distributions, Gurobi, Ipopt, RCall

import XLSX

function evsi(N,data1,data2,data3,line,demand,Safl)

    c = 0.75
    p_min = 2000
    p_max = 10000

    Xmax = 1000
    Dmax = demand*Xmax
    Y = 0.25*Safl*Xmax

    ben = 881

    alfa = zeros((N,N))
    alfa[1,:] = (ben/11)*[0 11 17 20]
    alfa[2,:] = (ben/11)*[11 0 22 33]
    alfa[3,:] = (ben/11)*[17 22 0 18]
    alfa[4,:] = (ben/11)*[20 33 18 0]

    valores = zeros((1,100))
    opt_gap = zeros((1,100))
    precios = zeros((N,100))
    concentraciones = zeros((N,100))
    flujos = zeros((N^2,100))
    dem_eff = zeros((N,100))
    infactible = zeros((1,100))

    x0 = data2[line,:]*Xmax
    d0 = data1[line,:]*Dmax

    t3 = time_ns()/1.0e9

        for sim=1:10

            m = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "OutputFlag"=>0, "TimeLimit" => 1800))
            @variable(m, v[1:N,1:N] >=0, Int)

            y = data3[sim,1:4]*Y
#            muest = data3[sim,5:8]

            delta = 0.9

            M = 4*p_max
            K = p_max*(3*Xmax + 1.3*Dmax)

            @variable(m, x[1:N] >=0)
            @variable(m, gamma[1:N])
            @variable(m, lambda[1:N,1:N])
            @variable(m, beta[1:N] >= 0)
            @variable(m, 0 <= w[1:N] <= 1, Int)
            @variable(m, 0 <= z[1:N,1:N] <= 1, Int)
            @variable(m, p_min <= p[1:N] <= p_max)
            @variable(m, 0 <= a[1:N] <= 1, Int)
            @variable(m, -K <= t1[1:N] <= K)
            @variable(m, -K <= t2[1:N] <= K)

            @expression(m, dEff[i=1:N], data3[sim,4+i]*d0[i]*(1-delta*((p[i]-p_min)/(p_max-p_min))))

            #Zonas: Santiago (1), Renca (2), Maipú (3), La Florida (4), Las Condes (5)
            #\alpha se calcula con distancia entre "centros" y tomando performance promedio: 11Km/Litro
            #Precio Promedio Bencina 95 en RM: $881 (20-05-2021)
            #Las distancias no son necesariamente simétricas, ¿considerar?

            @objective(m, Max, (1-c)*sum(p[i]*(-K*a[i] + t1[i]) + p[i]*dEff[i] for i=1:N))

            @constraint(m, [i=1:N],  t1[i] - K*a[i] <= 0)
            @constraint(m, [i=1:N],  t2[i] - K*(1 - a[i]) <= 0)
            @constraint(m, [i=1:N],  t1[i] >= 0)
            @constraint(m, [i=1:N],  t2[i] >= 0)

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
                    @constraint(m, c*(-beta[i]+beta[j]) - alfa[i,j] - lambda[i,j] + gamma[i] == 0)
                    #Desacople HC lambda
                    @constraint(m, 0 <= M*z[i,j] + lambda[i,j])
                    @constraint(m, 0 >= -M*z[i,j] + lambda[i,j])
                    #Desacople HC v
                    @constraint(m, 0 <= M*(1-z[i,j]) + v[i,j])
                    @constraint(m, 0 >= -M*(1-z[i,j]) + v[i,j])
                    #Desacople v-x0

                end
                #Desacople HC gammas
                @constraint(m, 0 <= gamma[i] + M*w[i])
                @constraint(m, 0 >= gamma[i] - M*w[i])

                @constraint(m, 0 <= sum(v[i,j] for j=1:N) - x0[i] + M*(1-w[i]))
                @constraint(m, sum(v[i,j] for j=1:N) - x0[i] - M*(1-w[i]) <= 0)

                @constraint(m, beta[i]-p[i] <= 0)

            end

            #Restricciones Cuadráticas
            @variable(m, 0 <= betaPrice[1:N] <= 1, Int)
            @variable(m, 0 <= betaCero[1:N] <= 1, Int)

            #If betaCero == 1 => beta== 0 y x-dEff >= 0
            @constraint(m,[i=1:N],  beta[i] <= p_max*(1-betaCero[i]) )
            @constraint(m,[i=1:N],  x[i]+y[i]-dEff[i] >= -K*(1-betaCero[i]) )

            #If betaPrice == 1 => beta== p y x-dEff <= 0
            @constraint(m,[i=1:N],  p[i]-beta[i] <= p_max*(1-betaPrice[i]) )
            @constraint(m,[i=1:N],  x[i]+y[i]-dEff[i] <= K*(1-betaPrice[i]) )

            #At most one side is active. If none is active, then x[i]-dEff[i] == 0
            @constraint(m, [i=1:N], betaPrice[i]+betaCero[i]<=1)
            @constraint(m, [i=1:N], x[i]+y[i] - dEff[i]<= K*(betaPrice[i]+betaCero[i]) );
            @constraint(m, [i=1:N], x[i]+y[i] - dEff[i]>= -K*(betaPrice[i]+betaCero[i]) );

            #######################################################

            status = optimize!(m)

            avance = string("EVSI", " Linea ","$line",", Dmax ", "$demand",", Y ","$Safl")

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

        cd("EV_Shared")

        CSV.write(string("Resultados_evsi1_",sufijo,".csv"), Data1)
        CSV.write(string("Resultados_evsi2_",sufijo,".csv"), Data2)
        CSV.write(string("Resultados_evsi3_",sufijo,".csv"), Data3)
        CSV.write(string("Resultados_evsi4_",sufijo,".csv"), Data4)
        CSV.write(string("Resultados_evsi5_",sufijo,".csv"), Data5)
        CSV.write(string("Resultados_evsi6_",sufijo,".csv"), Data6)

        cd("..")

end
