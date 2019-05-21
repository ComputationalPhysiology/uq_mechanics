#!/usr/bin/python




if __name__ == "__main__":
	from dolfin import *
	import chaospy as cp
        import numpy  as np
	import scipy 	 


	from ProblemSettings import *


	#random field
	from kle import KLE
	from covariance import CovarianceMatrix
	from kernel import Matern


	parameters['allow_extrapolation'] = True
	parameters["form_compiler"]["quadrature_degree"] = 3
	parameters["form_compiler"]["cpp_optimize"] = True
	parameters["form_compiler"]["representation"] = "uflacs"


	timer = Timer("00: forward run")
	model = problem()
        

	#b_ff = Constant(6.6)
	#b_xx = Constant(4.0)
	#b_fx = Constant(2.6)
	#CC   = Constant(1.1)      #kPa
	#K    = Constant(10.0)     #kPa

        #fiber_rotation_epi=50  # 50 
        #fiber_rotation_endo=40 # 40
        #sheet_rotation_epi=65  # 65
        #sheet_rotation_endo=25 # 25

	
        # Define Joint probability density funtion
        # joint = cp.J(alpha_dist,b_ff_dist,b_xx_dist,b_fx_dist,CC_dist,K_dist)

	# Number of terms in karhunen-loeve expansion
	#nkl = 65
	mu = 0.0
	sigma = 9.0
	eps1 = cp.Normal(mu = mu, sigma = sigma)
	eps2 = cp.Normal(mu = mu, sigma = sigma)
	eps3 = cp.Normal(mu = mu, sigma = sigma)
	eps4 = cp.Normal(mu = mu, sigma = sigma)
        eps5 = cp.Normal(mu = mu, sigma = sigma)
        eps6 = cp.Normal(mu = mu, sigma = sigma)
        eps7 = cp.Normal(mu = mu, sigma = sigma)
        eps8 = cp.Normal(mu = mu, sigma = sigma)
        eps9 = cp.Normal(mu = mu, sigma = sigma)
        eps10 = cp.Normal(mu = mu, sigma = sigma)
        eps11 = cp.Normal(mu = mu, sigma = sigma)
        eps12 = cp.Normal(mu = mu, sigma = sigma)
	eps13 = cp.Normal(mu = mu, sigma = sigma)
        eps14 = cp.Normal(mu = mu, sigma = sigma)
        eps15 = cp.Normal(mu = mu, sigma = sigma)
        eps16 = cp.Normal(mu = mu, sigma = sigma)
        eps17 = cp.Normal(mu = mu, sigma = sigma)
        eps18 = cp.Normal(mu = mu, sigma = sigma)
        eps19 = cp.Normal(mu = mu, sigma = sigma)
        eps20 = cp.Normal(mu = mu, sigma = sigma)
        eps21 = cp.Normal(mu = mu, sigma = sigma)
        eps22 = cp.Normal(mu = mu, sigma = sigma)
        eps23 = cp.Normal(mu = mu, sigma = sigma)
        eps24 = cp.Normal(mu = mu, sigma = sigma)
        eps25 = cp.Normal(mu = mu, sigma = sigma)
        eps26 = cp.Normal(mu = mu, sigma = sigma)
        eps27 = cp.Normal(mu = mu, sigma = sigma)
        eps28 = cp.Normal(mu = mu, sigma = sigma)
        eps29 = cp.Normal(mu = mu, sigma = sigma)
        eps30 = cp.Normal(mu = mu, sigma = sigma)
        eps31 = cp.Normal(mu = mu, sigma = sigma)
        eps32 = cp.Normal(mu = mu, sigma = sigma)
        eps33 = cp.Normal(mu = mu, sigma = sigma)
        eps34 = cp.Normal(mu = mu, sigma = sigma)
        eps35 = cp.Normal(mu = mu, sigma = sigma)
        eps36 = cp.Normal(mu = mu, sigma = sigma)
        eps37 = cp.Normal(mu = mu, sigma = sigma)
        eps38 = cp.Normal(mu = mu, sigma = sigma)
        eps39 = cp.Normal(mu = mu, sigma = sigma)
        eps40 = cp.Normal(mu = mu, sigma = sigma)
        eps41 = cp.Normal(mu = mu, sigma = sigma)
        eps42 = cp.Normal(mu = mu, sigma = sigma)
        eps43 = cp.Normal(mu = mu, sigma = sigma)
        eps44 = cp.Normal(mu = mu, sigma = sigma)
        eps45 = cp.Normal(mu = mu, sigma = sigma)
        eps46 = cp.Normal(mu = mu, sigma = sigma)
        eps47 = cp.Normal(mu = mu, sigma = sigma)
        eps48 = cp.Normal(mu = mu, sigma = sigma)
        eps49 = cp.Normal(mu = mu, sigma = sigma)
        eps50 = cp.Normal(mu = mu, sigma = sigma)
        eps51 = cp.Normal(mu = mu, sigma = sigma)
        eps52 = cp.Normal(mu = mu, sigma = sigma)
        eps53 = cp.Normal(mu = mu, sigma = sigma)
        eps54 = cp.Normal(mu = mu, sigma = sigma)
        eps55 = cp.Normal(mu = mu, sigma = sigma)
        eps56 = cp.Normal(mu = mu, sigma = sigma)
        eps57 = cp.Normal(mu = mu, sigma = sigma)
        eps58 = cp.Normal(mu = mu, sigma = sigma)
        eps59 = cp.Normal(mu = mu, sigma = sigma)
        eps60 = cp.Normal(mu = mu, sigma = sigma)
        eps61 = cp.Normal(mu = mu, sigma = sigma)
        eps62 = cp.Normal(mu = mu, sigma = sigma)
        eps63 = cp.Normal(mu = mu, sigma = sigma)
        eps64 = cp.Normal(mu = mu, sigma = sigma)
        eps65 = cp.Normal(mu = mu, sigma = sigma)

	joint_KL = cp.J(eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9,eps10,
			eps11,eps12,eps13,eps14,eps15,eps16,eps17,eps18,eps19,eps20,
			eps21,eps22,eps23,eps24,eps25,eps26,eps27,eps28,eps29,eps30,
			eps31,eps32,eps33,eps34,eps35,eps36,eps37,eps38,eps39,eps40,
			eps41,eps42,eps43,eps44,eps45,eps46,eps47,eps48,eps49,eps50,
			eps51,eps52,eps53,eps54,eps55,eps56,eps57,eps58,eps59,eps60,
			eps61,eps62,eps63,eps64,eps65)
	
        # Number of terms in karhunen-loeve expansion
        nkl = 65
	

        # Random Field
        p , l = np.inf, 1
        k = nkl
        kernel = Matern(p = p,l = l)    #0,1Exponential covariance kernel
        kle  = KLE(model.mesh, kernel, verbose = True)
        kle.compute_eigendecomposition(k = k)#20 by default
        #alpha_noise = kle.realizations(epsilon)
        #beta_noise = np.zeros(len(noise))


        # Monte Carlo	
        nsamples = 1
        qmc_scheme = "M"
        samples1= joint_KL.sample(nsamples, rule=qmc_scheme)
        
        
        pi = 0.0
        pf = 0.03
        dp = 0.01
        pressures = np.arange(pi, pf+dp, dp)

        #---------------------------------------------------
        print "---------------------------------"
        print "Quasi-Monte Carlo- Random Field"
        print "Sampling Method:", qmc_scheme
        print "Number of samples1:", nsamples
        print "Correlaton Length (cm):", l
	print "Sigma:", sigma
        print "Range of Pressures:", pressures
        print "---------------------------------"
        #---------------------------------------------------


        fileroot = "UQ_QMC"  + "_qmc" + str(qmc_scheme) + "_nmc"+ str(nsamples) + \
        "_pi" + str(int(pi*100)) + "_pf" + str(int(pf*100)) + "_dp" + str(int(dp*100)) +"_nKL" + str(nkl) + "_l" + str(l) + "_s" + str(int(sigma)) 

        qoi_names = np.array(["dVol", "dLen", "dRen", "dRep", "dThi", "dTwi"]).T

        txt_filename = str(fileroot) + "_qoi.txt"
        time_filename = str(fileroot) + "_time.txt"
        umean_filename = str(fileroot) + "_umean_.pvd"
        ustdv_filename  = str(fileroot) + "_ustd_.pvd"
        ums_filename = str(fileroot) + "_ums_.pvd"
        sffmean_filename = str(fileroot) + "_sffmean_.pvd"
        sffstdv_filename  = str(fileroot) + "_sffstd_.pvd"
	rfield_filename = str(fileroot) + "_rfield_.pvd"


        hist_filename = str(fileroot) + "_pressure" +  str(int(pf*100)) + "_hist.txt"


        txt_file = open(txt_filename,"w")
        time_file = open(time_filename, "w")
        hist_file = open(hist_filename,"w")
        umean_file = File(umean_filename)
        ustdv_file = File(ustdv_filename)
        ums_file = File(ums_filename)
        sffmean_file = File(sffmean_filename)
        sffstdv_file = File(sffstdv_filename)
	rfield_file = File(rfield_filename)


        umean = Function(model.V)
        ustdv = Function (model.V)
        ums = Function(model.V)
        sffmean = Function(model.projection_space)
        sffstdv = Function (model.projection_space)
        random_field = Function(FunctionSpace(model.mesh, "CG", 1))

        #model.update_params()
        i = 0
        qoi, u , sff = [], [], []
        for s in samples1.T:

          alpha_noise = kle.realizations(s)
          beta_noise = kle.realizations(s)
          model.update_randomfield( alpha_noise = alpha_noise, beta_noise = beta_noise )
          pressure = pi
          while pressure <= pressures[len(pressures)-1]:
                print "----------------"
                print "PRESSURE", pressure
                print "----------------"
                qoi_mc = model.run(pressure)

                pressure = pressure + dp


          qoi.append(qoi_mc[0:6])
          u.append(qoi_mc[6])
          sff.append(qoi_mc[7])
 
          i = i + 1
          if (i % nsamples/10 == 0):
              random_f = np.zeros((model.mesh.num_vertices()), dtype = float)
              random_f = alpha_noise
              random_field.vector()[:] = random_f[dof_to_vertex_map(FunctionSpace(model.mesh, "CG", 1))]
              rfield_file << random_field



        # output  qoi  metrics
        mean = np.mean(qoi, 0)
        std = np.std(qoi, 0)
        var = np.var(qoi, 0)
        cv = np.divide(np.array(std), np.array(mean))
        error = std / np.sqrt(nsamples)
        kurt = scipy.stats.kurtosis(qoi,0)
        skew = scipy.stats.skew(qoi,0)
        metrics = np.array([mean, std,var,cv,error,kurt,skew]).T
        rawdata = np.array(qoi[0:6])
	prodata =  np.column_stack((qoi_names, metrics.astype(np.object)))
	#prodata[:,1].astype(float)


        # saving qoi results
        txt_file.write("Pressure (kPa) = %f \n" % pressure)
        np.savetxt(txt_file, prodata, fmt=['%s','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f'], delimiter =  "  ",header = "qoi, mean, stdv, vari, cova, error, kurt, skew")
        #np.savetxt(txt_file, metrics, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f','%.5f'], header = "mean, stdv, vari, cova, error, kurt, skew")
        txt_file.write("\n")
        txt_file.write("\n")

        # saving qoi raw data
        np.savetxt(hist_file, rawdata, fmt=['%.10f','%.10f', '%.10f','%.10f','%.10f','%.10f'], header = "dVol, dLen, dRen, dRep, dThi, dTwi")


        # output  displacement metrics
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

        #closing files
        txt_file.close()
        time_file.close()
        hist_file.close()


