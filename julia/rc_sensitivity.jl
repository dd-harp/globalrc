using CSV
using Interpolations
using Random
using StatsPlots
using Turing


function infer_range(float_list)
    start = minimum(float_list)
    finish = maximum(float_list)
    count = length(float_list)
    delta = (finish - start) / (count - 1)
    start:delta:finish
end

infer_range(collect(0.1:.1:3))

for check in [0.01:.01:2.99, 1:1:5]
    @assert infer_range(collect(check)) == check
end


function interpolate_mesh()
    mesh_fn = "/home/adolgert/dev/analytics-pipeline/pr2ar_mesh.csv"
    row_cnt = 0
    for row in CSV.File(mesh_fn)
        row_cnt += 1
    end
    pr = zeros(Float64, row_cnt)
    rho = zeros(Float64, row_cnt)
    ar = zeros(Float64, row_cnt)
    row_idx = 1
    for row in CSV.File(mesh_fn)
        pr[row_idx] = row.PR
        rho[row_idx] = row.rho
        ar[row_idx] = row.AR
        row_idx += 1
    end
    rho_axis = sort(unique(rho))
    rho_scale = infer_range(rho_axis)
    pr_axis = sort(unique(pr))
    pr_scale = infer_range(pr_axis)
    rho_cnt = length(rho_axis)
    pr_cnt = length(pr_axis)
    @assert length(collect(pr_scale)) == pr_cnt
    @assert length(collect(rho_scale)) == rho_cnt
    rho_stride = findall(x -> x == rho_axis[1], rho)
    if rho_stride[1:10] == 1:10
        ar = reshape(ar, (rho_cnt, pr_cnt))
        println("rho rows, pr cols")
        itp = interpolate(ar, BSpline(Cubic(Line(OnGrid()))))
        sitp = Interpolations.scale(itp, rho_scale, pr_scale)
    else
        ar = reshape(ar, (pr_cnt, rho_cnt))
        itp = interpolate(ar, BSpline(Cubic(Line(OnGrid()))))
        println("pr rows, rho cols")
        sitp = Interpolations.scale(itp, pr_scale, rho_scale)
    end
    sitp
end

prfunc = interpolate_mesh()

# ## Turing model.
# @model function foi_model(x)
#     k ~ Beta(2, 5)
#     rho = k * x[1]
#     alpha = prfunc(x[2], rho)
#     x[3] = -log(1 - alpha) / 10
# end
# chn = sample(foi_model([0.3, 0.2, 0.01]), PG(10), 1000)


function foi_from_pram(pfpr, am, k = 0.6, tau = 10)
    rho = k * am
    alpha = prfunc(pfpr, rho)
    h = -log(1 - alpha) / tau
    h
end
foi_from_pram(0.4, 0.3)

@enum RcIn begin
    ampar=1  # multiplier on AM to get rho
    bpar  # biting
    cpar  #
    kpar  # The k in EIR
    rpar  # Recovery
    tpar  # days for foi
end
@enum RcOut alphap=1 kappap foip eirp vcp rcp

function pixel_work!(obs, pfpr, am, params)
    rho = params[Int(ampar)] * am
    alpha = prfunc(pfpr, rho)
    h = -log(1 - alpha) / params[Int(tpar)]
    kt = params[Int(tpar)] * params[Int(kpar)]
    eir = (exp(h * kt) - 1) / kt
    kappa = params[Int(cpar)] * pfpr
    V = eir / kappa
    D = params[Int(cpar)] / params[Int(rpar)]
    R = params[Int(bpar)] * V * D
    obs[Int(alphap)] = alpha
    obs[Int(kappap)] = kappa
    obs[Int(foip)] = h
    obs[Int(eirp)] = eir
    obs[Int(vcp)] = V
    obs[Int(rcp)] = R
end

params = zeros(Float64, Int(maximum(instances(RcIn))))
params[Int(ampar)] = 0.6
params[Int(bpar)] = 0.55
params[Int(cpar)] = .17
params[Int(kpar)] = 4.0
params[Int(rpar)] = 0.005
params[Int(tpar)] = 10
obs = zeros(Float64, Int(maximum(instances(RcOut))))
pixel_work!(obs, 0.3, 0.1, params)

function pin_pixel(pfpr, am)
    (obs, params) -> pixel_work!(obs, pfpr, am, params)
end
pin_pixel(0.3, 0.1)(obs, params)
using ForwardDiff
jacobian_buffer = ForwardDiff.JacobianConfig(pin_pixel(0.3, 0.1), obs, params)
fd = ForwardDiff.jacobian(pin_pixel(0.3, 0.1), obs, params, jacobian_buffer)
obs
H = fd' * fd
using LinearAlgebra
eigvals(H)
eigvecs(H)
fd