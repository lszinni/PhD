# ----------------------------------------------------------PACKAGES
using DelimitedFiles
using Printf
using TickTock
using Plots
using ColorSchemes
using LinearAlgebra

# ----------------------------------------------------------FUNCTIONS

function suma(ek, μ)
    sum = 0.0
    @inbounds @fastmath for i in eachindex(ek)
        sum += nf(ek[i] - μ, temp)
    end

    return 1.0 - δ - (2.0 * sum / length(ek))
    
end


function muCalculation(δ, t, tt, tpp, tz, temp, Nk, Nkz, Δ, ck, sk, cz, sz)
    
    Nd = 30
    μ = 0.0
    normak = 2.0 * Float64(Nkz) * (2.0 * Float64(Nk))^2    
    
    μA = -10.0
    μB = 10.0

    ek = [en(ck[ikx], ck[iky], cz[ikz], sk[ikx], sk[iky], sz[ikz], t, tt, tpp, tz, 0.0, δ, Δ) for ikx in eachindex(ck), iky in eachindex(ck), ikz in eachindex(cz)]    
    
    # suma(x) = 1.0 - δ - (2.0 * sum(nf.(ek .- x, temp)) / normak)
    
    sumaA = suma(ek, μA)
    sumaB = suma(ek, μB)
    
    for i in 1:Nd
        μ = (μA + μB) * 0.5

       sumaC = suma(ek, μ)

        if sumaA * sumaC > 0.0
            sumaA = sumaC
            μA = μ
        else
            sumaB = sumaC
            μB = μ
        end 
    end

    return μ
end

# ----------------------------------------------------------CALCULATE Δ

function deltaCalculation(t, tt, tz, tpp, temp, J, Nk, Nkz)
    
    normak = 2.0 * Float64(Nkz) * (2.0 * Float64(Nk))^2

    ck = cos.(Float64.(-Nk:Nk-1) * pi / Float64(Nk))
    sk = sin.(Float64.(-Nk:Nk-1) * pi / Float64(Nk))
    cz = cos.(Float64.(-Nkz:Nkz-1) * pi / Float64(Nkz))
    sz = sin.(Float64.(-Nkz:Nkz-1) * pi / Float64(Nkz))

    μ = 0.0
    xΔ = 0.0
    Δ = 1.0
    suma = 0.0

    while (abs(Δ - xΔ) > 1e-10)

        xΔ = Δ
        μ = muCalculation(δ, t, tt, tpp, tz, temp, Nk, Nkz, xΔ, ck, sk, cz, sz)
        
        ek = [en(ck[ikx], ck[iky], cz[ikz], sk[ikx], sk[iky], sz[ikz], t, tt, tpp, tz, μ, δ, xΔ) for ikx in eachindex(ck), iky in eachindex(ck), ikz in eachindex(cz)] 
        cksum = [ck[ikx] + ck[iky] for ikx in eachindex(ck), iky in eachindex(ck), ikz in eachindex(cz)] 
        # suma = sum(nf.(ek, temp) .* cksum)

        suma = 0.0
        @inbounds @fastmath for i in eachindex(ek)
            suma += nf(ek[i], temp) * cksum[i]
        end

        Δ = suma * J / (4.0 * normak) 
        
    end

    return Δ, μ

end

function nf(ek::Float64, temp::Float64)::Float64
    if temp > 0.0
        aa = ek / temp
        if aa > 60.0
            return 0.0
        elseif aa < -60.0
            return 1.0
        else
            return 1.0 / (exp(aa) + 1.0)
        end
    else
        return ek > 0.0 ? 0.0 : (ek < 0.0 ? 1.0 : 0.5)
    end
end

function en(cx, cy, cz, sx, sy, sz, t, tt, tpp, tz, μ, δ, Δ)
    fac1 = -(t * δ + 2.0 * Δ)
    fac2 = -2.0 * tt * δ
    # fac3 = -tz * δ
    fac3 = -0.25 * tz * δ * 0.5
    fac4 = -tpp * δ

    ekpar = fac1 * (cx + cy) + fac2 * cx * cy - μ
    ekperp = fac3 * ((cx - cy)^2) # * cz
    ekpp = fac4 * (cx^2 - sx^2 + cy^2 - sy^2)

    return ekpar + ekperp + ekpp 
end

function enmq(cx, cy, cz, sx, sy, sz, cxq, cyq, czq, sxq, syq, szq, t, tt, tpp, tz, μ, δ, Δ)
    fac1 = -(t * δ + 2.0 * Δ)
    fac2 = -2.0 * tt* δ
    # fac3 = -tz * δ
    fac3 = -0.25 * tz * δ * 0.5
    fac4 = -tpp * δ  

    ekpar = fac1 *  ((cx * cxq + sx * sxq) + (cy * cyq + sy * syq)) + fac2 * (cx * cxq + sx * sxq) * (cy * cyq + sy * syq) - μ
    ekperp = fac3 * (((cx * cxq + sx * sxq) - (cy * cyq + sy * syq))^2) # * (cz * czq + sz * szq)
    ekpp = fac4 * ((cx * cxq + sx * sxq)^2 - (sx * cxq - cx * sxq)^2 + (cy * cyq + sy * syq)^2 - (sy * cyq - cy * syq)^2)

    return ekpar + ekperp + ekpp 
end

function calcPi(t, tt, tz, tpp, temp, Γ, Δ, qx, qy, qz, wmax, μ, Nk , Nkz, Nw)

#  -----------Pi's and D order----------------
# 
#            R   L   rx  ry  Ax  Ay      
# 
#    R       1   2   3   4   5   6    
#    L           7   8   9   10  11   
#    rx              12  13  14  15
#    ry                  16  17  18
#    Ax                      19  20
#    Ay                          21 
# 
# --------------------------------------------
    normak = 2.0 * Float64(Nkz) * (2.0 * Float64(Nk))^2
    ω = Float64.(0:Nw) * wmax / Float64(Nw)
    pies = zeros(ComplexF64, length(ω), 3)

    ck = cos.(Float64.(-Nk:Nk-1) * pi / Float64(Nk))
    sk = sin.(Float64.(-Nk:Nk-1) * pi / Float64(Nk))
    cz = cos.(Float64.(-Nkz:Nkz-1) * pi / Float64(Nkz))
    sz = sin.(Float64.(-Nkz:Nkz-1) * pi / Float64(Nkz))
    cxq = cos(qx)
    cyq = cos(qy)
	czq = cos(qz)
    sxq = sin(qx)
    syq = sin(qy)
	szq = sin(qz)

    ek = [en(ck[ikx], ck[iky], cz[ikz], sk[ikx], sk[iky], sz[ikz], t, tt, tpp, tz, μ, δ, Δ) for ikx in eachindex(ck), iky in eachindex(ck), ikz in eachindex(cz)] 
    ekmq = [enmq(ck[ikx], ck[iky], cz[ikz], sk[ikx], sk[iky], sz[ikz], cxq, cyq, czq, sxq, syq, szq, t, tt, tpp, tz, μ, δ, Δ) for ikx in eachindex(ck), iky in eachindex(ck), ikz in eachindex(cz)] 
    ekk = [en(ck[ikx], ck[iky], cz[ikz], sk[ikx], sk[iky], sz[ikz], t, tt, tpp, tz, 0.0, δ, 0.0) for ikx in eachindex(ck), iky in eachindex(ck), ikz in eachindex(cz)] 
    ekkmq = [enmq(ck[ikx], ck[iky], cz[ikz], sk[ikx], sk[iky], sz[ikz], cxq, cyq, czq, sxq, syq, szq, t, tt, tpp, tz, 0.0, δ, 0.0) for ikx in eachindex(ck), iky in eachindex(ck), ikz in eachindex(cz)] 

    nfk = nf.(ek, temp)
    nfkmq = nf.(ekmq, temp)

    @inbounds @fastmath for i in eachindex(nfk)
        if abs(nfkmq[i] - nfk[i]) > 1e-13
            aux = 0.5 * (ekk[i] + ekkmq[i])
            denom = ω .+ im * Γ .- (ek[i] - ekmq[i])
            g = (nfkmq[i] .- nfk[i]) ./ denom
            pies[ : , 1] .+= g .* ((aux)^2) .- 0.5 .* (ekk[i] .- ekkmq[i]) .* nfk[i]
            pies[ : , 2] .+= g .* aux
            pies[ : , 3] .+= g 
        end
    end

    pies = -2.0 .* pies ./ normak

    return pies
end

function calcDab(t, tt, tz, tpp, temp, Γ, Δ, qx, qy, qz, wmax, μ, Nk , Nkz, Nw)

    Aq = α * (2.0 - cos(qx) - cos(qy)) + 1.0
    Vq = Vc / (Aq - cos(qz))
    Jq = 0.5 * J * (cos(qx) + cos(qy))
    @show Vq
    @show Jq
    println()

    print("Calculating Πab...")
    @time pies = calcPi(t, tt, tz, tpp, temp, Γ, Δ, qx, qy, qz, wmax, μ, Nk , Nkz, Nw)

    print("Calculating Dab...")
    ω = Float64.(0:Nw) * wmax / Float64(Nw)
    d = zeros(ComplexF64, length(ω), 3)
    dete = zeros(ComplexF64, length(ω))

   @time @inbounds @fastmath for i in eachindex(ω)
        dinv =  [
                δ^2 * (Vq - Jq) - pies[i, 1]    δ - pies[i, 2];
                δ - pies[i, 2]                  0.0 - pies[i, 3]
                ]
        dete[i] = det(dinv)
        dinv = inv(dinv)
        d[i, 1] = dinv[1, 1]
        d[i, 2] = dinv[1, 2]
        d[i, 3] = dinv[2, 2]        
    end

    pi11real = real.(pies[ : , 1])
    pi11im = imag.(pies[ : , 1])
    pi12real = real.(pies[ : , 2])
    pi12im = imag.(pies[ : , 2])
    pi22real = real.(pies[ : , 3])
    pi22im = imag.(pies[ : , 3])

    d11real = real.(d[ :, 1])
    d11im = imag.(d[ :, 1])
    d12real = real.(d[ :, 2])
    d12im = imag.(d[ :, 2])
    d22real = real.(d[ :, 3])
    d22im = imag.(d[ :, 3])

    detereal = real.(dete)
    deteim = imag.(dete)

    imax = argmax(d11im)

    return hcat(fill(qx, length(ω)), fill(qy, length(ω)), fill(qz, length(ω)), ω,
           pi11real, pi11im, pi12real, pi12im, pi22real, pi22im, 
           d11real, d11im, d12real, d12im, d22real, d22im,
           detereal, deteim),
           hcat(qx, qy, qz, ω[imax], pi11real[imax], pi11im[imax], pi12real[imax], pi12im[imax], pi22real[imax], pi22im[imax], 
           d11real[imax], d11im[imax], d12real[imax], d12im[imax], d22real[imax], d22im[imax],
           detereal[imax], deteim[imax])
end

function myprint(data, file_name)
    m = (a ->(@sprintf "%14.5e" a)).(data)
    writedlm(file_name, m, " ", quotes = false)
end


function main(qx, qy, qz, data_min, data_max)

    print("Calculating μ and Δ...")
    @time Δ, μ = deltaCalculation(t, tt, tz, tpp, temp, J, Nk, Nkz)
    @show μ
    @show Δ

    println()
    
    qx_actual = qx[data_min]
    qy_actual = qy[data_min]
    qz_actual = qz[data_min]

    # qx_actual = 0.1
    # qy_actual = 0.1
    # qz_actual = π

    println("\n---------------------------
        \nCalculating for       
        qx = $qx_actual,
        qy = $qy_actual,
        qz = $qz_actual,
        (data number: $data_min)    
    ")

    @time data, max_val = calcDab(t, tt, tz, tpp, temp, Γ, Δ, qx_actual, qy_actual, qz_actual, wmax, μ, Nk , Nkz, Nw)
    
    @time for i in data_min+1:data_max

        data_values = zeros(Float64, 1, 0)
        max_values = zeros(Float64, 1, 0)

        qx_actual = qx[i]
        qy_actual = qy[i]
        qz_actual = qz[i]

        println("\n---------------------------
            \nCalculating for       
            qx = $qx_actual,
            qy = $qy_actual,
            qz = $qz_actual,
            (data number: $i)    
        ")
            
        @time data_values, max_values = calcDab(t, tt, tz, tpp, temp, Γ, Δ, qx_actual, qy_actual, qz_actual, wmax, μ, Nk , Nkz, Nw)
    
        data = [data; data_values]
        max_val = [max_val; max_values]

    end

    myprint(data, "map")
    myprint(max_val, "max")

    # --------------PLOTING DATA

    x = qx[data_min:data_max]
    y = data[1:Nw+1, 4]
    z = data[:, 12]

    z_matrix = reshape(z, length(unique(y)), length(unique(x)))
    
    plot1 = heatmap(x, y, z_matrix, 
            c = :viridis, clim = (0, 50),
            xlabel = "qx",
            ylabel = "ω",
            xlim = (0, pi),
            ylim = (0, 2),
            )
    
    plot2 = plot!(x, max_val[:, 4],
        xlim = (0, pi),
        ylim = (0, 2),
        label = "max",
        color = "red")

    fig = scatter!(x, max_val[:, 4],
        xlim = (0, pi),
        ylim = (0, 2),
        markersize = 2,
        color = "red",
        primary = false)

    savefig(fig, "fig_largeN_test.pdf")
end




# ----------------------------------------------------------PARAMETERS

Nk = 400
Nkz = 1
Nw = 4000

t = 0.8
tt = -0.09 * t 
tpp = 0.07 * t
tz = 0.01 * t
J = 0.3 * t
δ = 0.16
temp = 0.0001
Γ = 0.04
Vc = 18.8 
α = 4.1 
wmax = 10.0 * t
# meff = 1.0

# ----------------------------------------------------------PRINT INFO

println("Parameters for calculation:
        δ = $δ
        t = $t
        t'= $tt
        t'' = $tpp
        tz = $tz
        J = $J
        Temp = $temp
        Γ = $Γ
        Vc = $Vc
        α = $alpha

        Nk = $Nk
        Nkz = $Nkz
        Nw = $Nw
")

# ----------------------------------------------------------DATA GENERATION

Nqx = 50
Nqy = 1
Nqz = 1

qx = Float64.(0:Nqx) * pi / Float64(Nqx)
# qy = Float64.(0:Nqy) * pi / Float64(Nqy)
# qz = Float64.(0:Nqz) * pi / Float64(Nqz)

qy = zeros(Float64, length(qx))
qz = fill(pi, length(qx))

# ----------------------------------------------------------SET MIN AND MAX VALUES OF DATA

data_min = 1

data_max = 51

# ----------------------------------------------------------START CALCULATION

println("Calculation from $data_min to $data_max")

tick()

main(qx, qy, qz, data_min, data_max)

tock()
