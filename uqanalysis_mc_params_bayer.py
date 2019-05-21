#!/usr/bin/python


from dolfin import *
import numpy  as np
import scipy  
import chaospy as cp
from ProblemSettings import *


parameters['allow_extrapolation'] = True
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"



timer = Timer("00: forward run")
model = problem()
        
# Define random distribution for each parameter
b_ff_dist = cp.Normal(mu = 6.6, sigma = 0.99)
b_xx_dist = cp.Normal(mu = 4.0, sigma = 0.6)
b_fx_dist = cp.Normal(mu = 2.6, sigma = 0.39)
CC_dist   = cp.Lognormal (mu = -0.049, sigma = 0.538)#cp.Normal(mu = 1.1 , sigma = 0.165)      #kPa
K_dist    = cp.Lognormal(mu = 2.2964, sigma = 0.11)#cp.Normal( mu = 10.0, sigma = 1.5)           #kPa
fiber_rotation_epi_dist = cp.Normal(mu = 50,sigma = 7.5)                 #deg
fiber_rotation_end_dist = cp.Normal(mu = 40,sigma = 6.0)                 #deg
sheet_rotation_epi_dist = cp.Normal(mu = 65,sigma = 9.75)                #deg
sheet_rotation_end_dist = cp.Normal(mu = 25,sigma = 3.75)                #deg

print('Calling Chaospy')

# Define Joint probability density funtion
#joint = cp.J(b_ff_dist,b_xx_dist,b_fx_dist,CC_dist,K_dist,fiber_rotation_epi_dist)
joint = cp.J(b_ff_dist,b_xx_dist,b_fx_dist,CC_dist,K_dist,fiber_rotation_epi_dist,fiber_rotation_end_dist,sheet_rotation_epi_dist,sheet_rotation_end_dist)


print(joint)

print('done with Chaospy')

# Monte Carlo	
nsamples = 5
qmc_scheme = "H"
samples= joint.sample(nsamples, rule=qmc_scheme)
	
# range of pressure
pi = 0.0
pf = 2.0
dp = 0.05
pressures = np.arange(pi, pf+dp, dp)


        
# --------------------------------------------------
print("---------------------------------")
print("Quasi-Monte Carlo: parameters")
print("Sampling Method:", qmc_scheme)
print("Number of samples:", nsamples)
print("Range of Pressures:", pressures)
print("---------------------------------")
#---------------------------------------------------


# open files
fileroot = "UQ_QMC"  + "_qmc" + str(qmc_scheme) + "_nmc"+ str(nsamples) + \
  "_pi" + str(int(pi*100)) + "_pf" + str(int(pf*100)) + "_dp" + str(int(dp*100)) +\
  "_sp" + str(int(len(joint))) +"_allqoi"

qoi_names = np.array(["dVol"])#, "dLen","dLep","dL", "dRen", "dRep", "dThi", "dTwi","dWV"]).T


txt_filename = str(fileroot) + "_qoi.txt"
time_filename = str(fileroot) + "_time.txt" 
umean_filename = str(fileroot) + "_umean_.pvd"
ustdv_filename  = str(fileroot) + "_ustd_.pvd"
ums_filename = str(fileroot) + "_ms_.pvd"
sffmean_filename  = str(fileroot) + "_sffmean_.pvd"
sffstdv_filename  = str(fileroot) + "_sffstdv_.pvd"


ums_filename = str(fileroot) + "_ms_.pvd"
hist_filename = str(fileroot) + "_pressure" +  str(int(pf*100)) + "_hist.txt"


txt_file = open(txt_filename,"w")
time_file = open(time_filename, "w")
hist_file = open(hist_filename,"w")

umean_file = File(umean_filename)
ustdv_file = File(ustdv_filename)
ums_file = File(ums_filename)
sffmean_file =  File(sffmean_filename)
sffstdv_file = File(sffstdv_filename)


umean = Function(model.V)
ustdv = Function (model.V)
ums = Function(model.V)
sffmean = Function(model.projection_space) 
sffstdv = Function(model.projection_space)
	
nnan,i = 0,0
qoi, u , sff = [], [], []
#model = problem()#WARNING!!
for s in samples.T:
    model = problem()
    model.update_params(*s)
    pressure = pi
    while pressure <= pressures[len(pressures)-1]:
        print("----------------")
        print("PRESSURE", pressure)
        print("----------------")

        qoi_mc = model.run(pressure)
        if (i == 0) and (pressure == 0.0): 
            txt_file.write("--------------------------------- \n")
            txt_file.write("QMC: parameters \n")
            txt_file.write("Number of samples: %i \n" % nsamples)
            txt_file.write("Sampling Method: %s \n" % qmc_scheme)
            txt_file.write("Range of Pressures: %.2f,%.2f,%.2f \n" % (pi ,pf ,dp ))
            txt_file.write("--------------------------------- \n")
            txt_file.write("\n")
            txt_file.write("Initial Values: \n" )
            #np.savetxt(txt_file,  np.column_stack((qoi_names, np.array(qoi_mc[0:9]).T.astype(np.object))) , fmt=['%s','%.8f'], delimiter =  "  ")
            txt_file.write("\n")
            txt_file.write("\n")

        while np.isnan(qoi_mc[0]) == True:
            print("WARNING: NaN")
            newsample = joint.sample(1, rule = "R")
            model = problem()
            model.update_params(*newsample.T[0])
            nnan = nnan + 1

        #pressure = pi
        qoi_mc = model.run(pressure)

        pressure = pressure + dp
          
    qoi.append(qoi_mc[0:9])
    u.append(qoi_mc[9])
    sff.append(qoi_mc[10])
    i = i +1
    print("-------------------------")
    print("Number of MC steps:", i)
    print("-------------------------")
          
# output  qoi  metrics
mean = np.mean(qoi, 0)
std = np.std(qoi, 0)
var = np.var(qoi, 0)
cv = np.divide(np.array(std), np.array(mean))
error = std / np.sqrt(nsamples)
kurt = scipy.stats.kurtosis(qoi,0)
skew = scipy.stats.skew(qoi,0)
metrics = np.array([mean, std,var,cv,error,kurt,skew]).T
rawdata = np.array(qoi)
prodata =  np.column_stack((qoi_names, metrics.astype(np.object)))
#prodata[:,1].astype(float)

# saving qoi results
pressure = pressure - dp
txt_file.write("Pressure (kPa) = %f \n" % pressure)
np.savetxt(txt_file, prodata, fmt=['%s','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f'], delimiter =  "  ",header = "qoi, mean, stdv, vari, cova, error, kurt, skew")
txt_file.write("\n")
txt_file.write("\n")
#np.savetxt(txt_file, metrics, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f','%.5f'], header = "mean, stdv, vari, cova, error, kurt, skew")

# saving qoi raw data
np.savetxt(hist_file, rawdata, fmt=['%.10f','%.10f','%.10f','%.10f' ,'%.10f','%.10f','%.10f','%.10f','%.10f'], \
               header = "dVol, dLen, dLep, dL, dRen, dRep, dThi, dTwi, dWV")

# output  displacement/stress metrics
u_mean = np.mean(u, 0)
u_stdv = np.std(u, 0)
sff_mean = np.mean(sff,0)
sff_stdv = np.std(sff,0)

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
time_file.write("Time (sec) %s \n" %time)
time_file.write("\n")
time_file.write("Number of NaN: %i" %nnan)

#closing files
txt_file.close()
time_file.close()
hist_file.close()




	
	
