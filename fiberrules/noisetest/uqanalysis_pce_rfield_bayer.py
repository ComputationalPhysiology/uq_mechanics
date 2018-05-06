#!/usr/bin/python




if __name__ == "__main__":
	from dolfin import *
	import chaospy as cp
        import numpy  as np

	from ProblemSettings import *


        #random field
        from kle import KLE
        from covariance import CovarianceMatrix
        from kernel import Matern



	parameters["form_compiler"]["quadrature_degree"] = 3
	parameters["form_compiler"]["cpp_optimize"] = True
	parameters["form_compiler"]["representation"] = "uflacs"


	timer = Timer("00: forward run")
	model = problem()



	# Parameters are 
	#b_ff = Constant(6.6)
	#b_xx = Constant(4.0)
	#b_fx = Constant(2.6)
	#CC   = Constant(1.1)      #kPa
	#K    = Constant(10.0)     #kPa
	#alpha = 60.0


	
        # Number of terms in karhunen-loeve expansion
        #nkl = 10
	mu = 0.0
	sigma = 1.0  
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

        joint_KL = cp.J(eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9,eps10)
                       # eps11,eps12,eps13,eps14,eps15,eps16,eps17,eps18,eps19,eps20,
                       # eps21,eps22,eps23,eps24,eps25,eps26,eps27,eps28,eps29,eps30,
                       # eps31,eps32,eps33,eps34,eps35,eps36,eps37,eps38,eps39,eps40)
                       # eps41,eps42,eps43,eps44,eps45,eps46,eps47,eps48,eps49,eps50,
                       # eps51,eps52,eps53,eps54,eps55,eps56,eps57,eps58,eps59,eps60,
                       # eps61,eps62,eps63,eps64,eps65)

        # Number of terms in karhunen-loeve expansion
        nkl = 10

        # Polynomial Chaos Expansion: three terms recursion relation
        order = 1
        Polynomials = cp.orth_ttr(order, joint_KL)

        # Generate samples by quasi-random samples
        qmc_scheme = "H"
        nnodes = 2*len(Polynomials)
        nodes = joint_KL.sample(nnodes, qmc_scheme)

        # Regression method in Point Collocation
        rgm_rule ="T" #Orthogonal Matching Pursuit#"LS"


        # Random Field
        p , l = np.inf, 1
        k = nkl
        kernel = Matern(p = p,l = l)    #0,1Exponential covariance kernel
        kle  = KLE(model.mesh, kernel, verbose = True)
        kle.compute_eigendecomposition(k = k)#20 by default
        #real = kle.realizations(epsilon)

	# range of pressure
        pi = 0.0
        pf = 0.03
        dp = 0.01
        pressures = np.arange(pi, pf+dp, dp)

        print "---------------------------------"
        print "PCE-Random Field"
        print "Nkl terms:", int(nkl)
        print "Order:", order
        print "Sampling Method:", qmc_scheme
        print "Number of nodes:", nnodes
        print "Regression method:", rgm_rule
        print "Correlaton Length:", l
        print "Sigma:", sigma
        print "Range of Pressures:", pressures
        print "---------------------------------"



        # open files
        fileroot = "UQ_PCE-PC"  + "_o" + str(order) + "_qmc" + str(qmc_scheme) + "_nnodes"+ str(nnodes) +  "_rgm" + str(rgm_rule) +\
        "_pi" + str(int(pi*100)) + "_pf" + str(int(pf*100)) + "_dp" + str(int(dp*100)) + "_nKL" + str(nkl) + "_l" + str(l) + "_s" + str(int(sigma))

        qoi_names = np.array(["dVol", "dLen", "dRen", "dRep", "dThi", "dTwi"]).T

        txt_filename = str(fileroot) + "_qoi.txt"
        time_filename = str(fileroot) + "_time.txt"
        umean_filename = str(fileroot) + "_mean_.pvd"
        ustdv_filename  = str(fileroot) + "_std_.pvd"
        ums_filename = str(fileroot) + "_ms_.pvd"
        rfield_filename = str(fileroot) + "_rfield_.pvd"
        hist_filename = str(fileroot) +   "_hist.txt"
        sffmean_filename = str(fileroot) + "_sffmean.pvd"
        sffstdv_filename  = str(fileroot) + "_sffstd_.pvd"

        sai_filename = str(fileroot) + "_SAi_.txt"
        usai1_filename =str(fileroot) + "_pressure" +  str(int(pf*100)) + "_uSAi1_.pvd"
        usait_filename =str(fileroot) + "_pressure" +  str(int(pf*100)) + "_uSAit_.pvd"
        sffsai1_filename =str(fileroot) + "_pressure" +  str(int(pf*100)) + "_sffSAi1_.pvd"
        sffsait_filename =str(fileroot) + "_pressure" +  str(int(pf*100)) + "_sffSAit_.pvd"
        hist_filename = str(fileroot) + "_pressure" +  str(int(pf*100)) + "_hist.txt"




        # Open file
        hist_file = open(hist_filename,"w")
        txt_file = open(txt_filename,"w")
        sai_file = open(sai_filename,"w")
        time_file = open(time_filename, "w")


        umean_file = File(umean_filename)
        ustdv_file = File(ustdv_filename)
        ums_file = File(ums_filename)
        usai1_file = File(usai1_filename)
        usait_file = File(usait_filename)
        rfield_file = File(rfield_filename)

        sffmean_file = File(sffmean_filename)
        sffstdv_file = File(sffstdv_filename)
        sffsai1_file = File(sffsai1_filename)
        sffsait_file = File(sffsait_filename)




        umean = Function(model.V)
        ustdv = Function (model.V)
        ums = Function(model.V)
        usai1 = Function(model.V)
        usait = Function(model.V)

        sffmean = Function(model.projection_space)
        sffstdv = Function(model.projection_space)
        sffsai1 = Function(model.projection_space)
        sffsait = Function(model.projection_space)

        random_field = Function(FunctionSpace(model.mesh, "CG", 1))


        # Evaluate model in sample points
        solves1 = []
        solves2 = []
        solves3 = []
        nnan = 0
        for s in nodes.T:
          alpha_noise = kle.realizations(s)
          model.update_randomfield( alpha_noise = alpha_noise)

          pressure = pi
          while pressure <= pressures[len(pressures)-1]:
                print "----------------"
                print "PRESSURE", pressure
                print "----------------"
                qoi = model.run(pressure)

                while np.isnan(qoi[0]) == True:
                        print "WARNING: NaN"
                        newsample = joint.sample(1, rule = "R")
                        model = problem()
                        model.update_params(*newsample.T[0])
                        nnan = nnan + 1
                        pressure = pi

                pressure = pressure + dp


          solves1.append(qoi[0:6])
          solves2.append(qoi[6])
          solves3.append(qoi[7])








        # Create surrogate model
        #u_hat(x,t;q)=sum_n C_n(x) * P_n(q)
        qoi1_hat =  cp.fit_regression(Polynomials, nodes, solves1,rule=rgm_rule)
        qoi2_hat = cp.fit_regression(Polynomials, nodes, solves2,rule=rgm_rule)
        qoi3_hat =  cp.fit_regression(Polynomials, nodes, solves3,rule=rgm_rule)

       #output qoi metrics
        mean = cp.E(qoi1_hat, joint_KL)
        var = cp.Var(qoi1_hat, joint_KL)
        std = np.sqrt(var)
        cv = np.divide(np.array(std), np.array(mean))
        kurt = cp.Kurt(qoi1_hat,joint_KL)
        skew = cp.Skew(qoi1_hat,joint_KL)
        metrics = np.array([mean, std,var,cv,kurt,skew]).T
        rawdata = np.array(qoi1_hat(*joint_KL.sample(10**5))).T#np.array(solves1)
        prodata =  np.column_stack((qoi_names, metrics.astype(np.object)))
        #prodata[:,1].astype(float)


        # saving qoi results
        txt_file.write("Pressure (kPa) = %f \n" % pressure)
        np.savetxt(txt_file, prodata, fmt=['%s','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f'], delimiter =  "  ",header = "qoi, mean, stdv, vari, cova, kurt, skew")
        txt_file.write("\n")
        txt_file.write("\n")
        #np.savetxt(txt_file, metrics,fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f'], header = "mean, stdv, vari, cova, kurt, skew")
        # saving qoi raw data
        np.savetxt(hist_file, rawdata, fmt=['%.10f','%.10f', '%.10f','%.10f','%.10f','%.10f'], delimiter = " ", header = "dVol, dLen, dRen, dRep, dThi, dTwi")

        # output sensitivity index 
        sai1 = cp.Sens_m(qoi1_hat,joint_KL)
        sait = cp.Sens_t(qoi1_hat,joint_KL)
        # save sensitivity index 
        sai_file.write("Pressure (kPa) = %f \n" % pressure)
        sai_file.write("First order Sensitivity index \n" )
        #np.savetxt(sai_file, sai1.T, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f'], header = " b_ff, b_xx, b_fx, C, K, fr_epi")
        #np.savetxt(sai_file, sai1.T, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f','%.5f','%.5f','%.5f'], header = " b_ff, b_xx, b_fx, C, K, fr_epi, fr_endo, sr_epi, sr_endo")

        sai_file.write("\n")
        sai_file.write("Total-effect index \n" )
        #np.savetxt(sai_file, sait.T, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f'], header = " b_ff, b_xx, b_fx, C, K, fr_epi")
        #np.savetxt(sai_file, sait.T, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f','%.5f','%.5f','%.5f'], header = " b_ff, b_xx, b_fx, C, K, fr_epi,fr_endo, sr_epi, sr_endo")

        sai_file.write("\n")
        sai_file.write("\n")

        #output displacement/stress metrics
        #sai1 = cp.Sens_m(qoi_hat[0],joint)
        u_mean = cp.E(qoi2_hat, joint_KL)
        u_var = cp.Var(qoi2_hat, joint_KL)
        u_stdv = np.sqrt(u_var)
        u_sai1 = cp.Sens_m(qoi2_hat,joint_KL)
        u_sai2 = cp.Sens_m2(qoi2_hat,joint_KL)
        u_sait = cp.Sens_t(qoi2_hat,joint_KL)

        sff_mean = cp.E(qoi3_hat, joint_KL)
        sff_var = cp.Var(qoi3_hat, joint_KL)
        sff_stdv = np.sqrt(sff_var)
        sff_sai1 = cp.Sens_m(qoi3_hat,joint_KL)
        sff_sai2 = cp.Sens_m2(qoi3_hat,joint_KL)
        sff_sait = cp.Sens_t(qoi3_hat,joint_KL)


        # save displacement/stress results
        umean.vector()[:] = u_mean
        ustdv.vector()[:] = u_stdv
        ums.vector()[:] = u_mean + u_stdv
        umean_file << umean
        ustdv_file << ustdv
        ums_file << ums

        sffmean.vector()[:] = sff_mean
        sffstdv.vector()[:] = sff_stdv
        sffmean_file << sffmean
        sffstdv_file << sffstdv


        for param in range(len(nodes)):
                usai1.vector()[:] = u_sai1[param]
                usait.vector()[:] = u_sait[param]
                usai1_file << usai1
                usait_file << usait
                sffsai1.vector()[:] = sff_sai1[param]
                sffsait.vector()[:] = sff_sait[param]
                sffsai1_file << sffsai1
                sffsait_file << sffsait



          #else: continue



        time = timer.stop()
        list_timings(TimingClear_clear, [TimingType_wall, TimingType_user, TimingType_system])
        # save total user time
        time_file.write("Time (sec) %s \n" %time)
        time_file.write("\n")
        time_file.write("Number of NaN: %i" %nnan)


        #closing files
        txt_file.close()
        time_file.close()
        sai_file.close()
        hist_file.close()
