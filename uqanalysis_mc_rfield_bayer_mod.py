from dolfin import *
import chaospy as cp
import numpy as np
import scipy

from ProblemSettings import *

# random field
from kle import KLE
from covariance import CovarianceMatrix
from kernel import Matern

parameters['allow_extrapolation'] = True
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"

timer = Timer("00: forward run")
model = problem()

# b_ff = Constant(6.6)
# b_xx = Constant(4.0)
# b_fx = Constant(2.6)
# CC   = Constant(1.1)      #kPa
# K    = Constant(10.0)     #kPa

# fiber_rotation_epi=50  # 50
# fiber_rotation_endo=40 # 40
# sheet_rotation_epi=65  # 65
# sheet_rotation_endo=25 # 25


# Define Joint probability density funtion
# joint = cp.J(alpha_dist,b_ff_dist,b_xx_dist,b_fx_dist,CC_dist,K_dist)

# Number of terms in karhunen-loeve expansion
# nkl = 65
mu = 0.0
s = 1.0
eps1 = cp.Normal(mu=mu, sigma=s)
eps2 = cp.Normal(mu=mu, sigma=s)
eps3 = cp.Normal(mu=mu, sigma=s)
eps4 = cp.Normal(mu=mu, sigma=s)
eps5 = cp.Normal(mu=mu, sigma=s)
eps6 = cp.Normal(mu=mu, sigma=s)
eps7 = cp.Normal(mu=mu, sigma=s)
eps8 = cp.Normal(mu=mu, sigma=s)
eps9 = cp.Normal(mu=mu, sigma=s)
eps10 = cp.Normal(mu=mu, sigma=s)
eps11 = cp.Normal(mu=mu, sigma=s)
eps12 = cp.Normal(mu=mu, sigma=s)
eps13 = cp.Normal(mu=mu, sigma=s)
eps14 = cp.Normal(mu=mu, sigma=s)
eps15 = cp.Normal(mu=mu, sigma=s)
eps16 = cp.Normal(mu=mu, sigma=s)
eps17 = cp.Normal(mu=mu, sigma=s)
eps18 = cp.Normal(mu=mu, sigma=s)
eps19 = cp.Normal(mu=mu, sigma=s)
eps20 = cp.Normal(mu=mu, sigma=s)

joint_KL = cp.J(eps1, eps2, eps3, eps4, eps5, eps6, eps7, eps8, eps9, eps10,
                eps11, eps12, eps13, eps14, eps15, eps16)  # ,eps17,eps18,eps19,eps20)
# eps21,eps22,eps23,eps24,eps25,eps26,eps27,eps28,eps29,eps30,
# eps31,eps32,eps33,eps34,eps35,eps36,eps37,eps38,eps39,eps40,
# eps41,eps42,eps43,eps44,eps45,eps46,eps47,eps48,eps49,eps50,
# eps51,eps52,eps53,eps54,eps55,eps56,eps57,eps58,eps59,eps60,
# eps61,eps62,eps63,eps64,eps65)

# Number of terms in karhunen-loeve expansion
nkl = 16

# Random Field
p, l = np.inf, 3
k = nkl
sigma = 10 * np.pi / 180.
kernel = Matern(p=p, l=l, s=sigma)  # 0,1Exponential covariance kernel
kle = KLE(model.mesh, kernel, verbose=True)
kle.compute_eigendecomposition(k=k)  # 20 by default
# alpha_noise = kle.realizations(epsilon)
# beta_noise = np.zeros(len(noise))


# Monte Carlo
nsamples = 30
qmc_scheme = "H"
samples1 = joint_KL.sample(nsamples, rule=qmc_scheme)

pi = 0.0
pf = 2.0
dp = 0.05
pressures = np.arange(pi, pf + dp, dp)


print("---------------------------------")
print("Quasi-Monte Carlo- Random Field")
print("Sampling Method:", qmc_scheme)
print("Number of samples:", nsamples)
print("Correlaton Length (cm):", l)
print("Sigma:", sigma)
print("Range of Pressures:", pressures)
print("---------------------------------")

fileroot = "UQ_QMC" + "_qmc" + str(qmc_scheme) + "_nmc" + str(nsamples) + \
           "_pi" + str(int(pi * 100)) + "_pf" + str(int(pf * 100)) + "_dp" + str(int(dp * 100)) + \
           "_nKL" + str(nkl) + "_l" + str(l) + "_s" + str(int(sigma)) + \
           "_srf" + "_new"

qoi_names = np.array(["dVol", "dLen", "dLep", "dL", "dRen", "dRep", "dThi", "dTwi", "dWV"]).T

txt_filename = str(fileroot) + "_qoi.txt"
time_filename = str(fileroot) + "_time.txt"
umean_filename = str(fileroot) + "_umean_.pvd"
ustdv_filename = str(fileroot) + "_ustd_.pvd"
ums_filename = str(fileroot) + "_ums_.pvd"
sffmean_filename = str(fileroot) + "_sffmean_.pvd"
sffstdv_filename = str(fileroot) + "_sffstd_.pvd"
rfield_filename = str(fileroot) + "_rfield_.pvd"

hist_filename = str(fileroot) + "_pressure" + str(int(pf * 100)) + "_hist.txt"

txt_file = open(txt_filename, "w")
time_file = open(time_filename, "w")
hist_file = open(hist_filename, "w")
umean_file = File(umean_filename)
ustdv_file = File(ustdv_filename)
ums_file = File(ums_filename)
sffmean_file = File(sffmean_filename)
sffstdv_file = File(sffstdv_filename)
rfield_file = File(rfield_filename)

umean = Function(model.V)
ustdv = Function(model.V)
ums = Function(model.V)
sffmean = Function(model.projection_space)
sffstdv = Function(model.projection_space)
random_field = Function(FunctionSpace(model.mesh, "CG", 1))

# model.update_params()
i, nnan = 0, 0
qoi, u, sff = [], [], []
for s in samples1.T:

    alpha_noise = kle.realizations(s)[dof_to_vertex_map(model.fiber_space)]
    beta_noise = 0.0  # kle.realizations(s)
    model.update_randomfield(alpha_noise=alpha_noise * 180. / np.pi, beta_noise=beta_noise * 180. / np.pi)
    pressure = pi
    while pressure <= pressures[len(pressures) - 1]:
        print("----------------")
        print("PRESSURE", pressure)
        print("----------------")
        qoi_mc = model.run(pressure)

        if (i == 0) & (pressure == 0.0):
            txt_file.write("--------------------------------- \n")
            txt_file.write("QMC: Random Field \n")
            txt_file.write("Nkl terms: %i \n" % int(nkl))
            txt_file.write("Number of samples: %i \n" % nsamples)
            txt_file.write("Sampling Method: %s \n" % qmc_scheme)
            txt_file.write("Correlaton Length: %i" % l)
            txt_file.write("Sigma: %i \n" % int(sigma))
            txt_file.write("Range of Pressures: %.2f,%.2f,%.2f \n" % (pi, pf, dp))
            txt_file.write("--------------------------------- \n")
            txt_file.write("\n")

            txt_file.write("Initial Values: \n")
            np.savetxt(txt_file, np.column_stack((qoi_names, np.array(qoi_mc[0:9]).T.astype(np.object))),
                       fmt=['%s', '%.8f'], delimiter="  ")
            txt_file.write("\n")
            txt_file.write("\n")

        while np.isnan(qoi_mc[0]) == True:
            print("WARNING: NaN")
            newsample = joint_KL.sample(1, rule="R")
            model = problem()
            alpha_noise = kle.realizations(newsample.T[0])[dof_to_vertex_map(model.fiber_space)]
            model.update_randomfield(alpha_noise=alpha_noise)
            nnan = nnan + 1
            pressure = pi
            qoi_mc = model.run(pressure)

        pressure = pressure + dp

    qoi.append(qoi_mc[0:9])
    u.append(qoi_mc[9])
    sff.append(qoi_mc[10])

    i = i + 1
    print("Number of MC steps:", i)
    if (i % nsamples / 10 == 0):
        random_f = np.zeros((model.mesh.num_vertices()), dtype=float)
        random_f = alpha_noise
        random_field.vector()[:] = random_f[dof_to_vertex_map(FunctionSpace(model.mesh, "CG", 1))]
        rfield_file << random_field

# output  qoi  metrics
mean = np.mean(qoi, 0)
std = np.std(qoi, 0)
var = np.var(qoi, 0)
cv = np.divide(np.array(std), np.array(mean))
error = std / np.sqrt(nsamples)
kurt = scipy.stats.kurtosis(qoi, 0)
skew = scipy.stats.skew(qoi, 0)
metrics = np.array([mean, std, var, cv, error, kurt, skew]).T
rawdata = np.array(qoi[0:9])
prodata = np.column_stack((qoi_names, metrics.astype(np.object)))
# prodata[:,1].astype(float)


# saving qoi results
txt_file.write("Pressure (kPa) = %f \n" % pressure)
np.savetxt(txt_file, prodata, fmt=['%s', '%.8f', '%.8f', '%.8f', '%.8f', '%.8f', '%.8f', '%.8f'], delimiter="  ",
           header="qoi, mean, stdv, vari, cova, error, kurt, skew")
# np.savetxt(txt_file, metrics, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f','%.5f'], header = "mean, stdv, vari, cova, error, kurt, skew")
txt_file.write("\n")
txt_file.write("\n")

# saving qoi raw data
np.savetxt(hist_file, rawdata, fmt=['%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f'], \
           header="dVol, dLen, dLep, dL, dRen, dRep, dThi, dTwi, dWV")

# output  displacement metrics
u_mean = np.mean(u, 0)
u_stdv = np.std(u, 0)
sff_mean = np.mean(sff, 0)
sff_stdv = np.std(sff, 0)

# save displacement results
umean.vector()[:] = u_mean
ustdv.vector()[:] = u_stdv
ums.vector()[:] = u_mean + u_stdv
umean_file << umean
ustdv_file << ustdv
ums_file << ums
# save stress projection
sffmean.vector()[:] = sff_mean
sffstdv.vector()[:] = sff_stdv
sffmean_file << sffmean
sffstdv_file << sffstdv

# stop time
time = timer.stop()
list_timings(TimingClear_clear, [TimingType_wall, TimingType_user, TimingType_system])

# save total user time
time_file.write("Time (sec) %s \n" % time)

# closing files
txt_file.close()
time_file.close()
hist_file.close()
