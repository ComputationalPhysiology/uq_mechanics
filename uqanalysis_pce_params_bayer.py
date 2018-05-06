#!/usr/bin/python




if __name__ == "__main__":
	from dolfin import *
        import numpy  as np
	import chaospy as cp
        import pandas as pd

	from ProblemSettings import *



	parameters['allow_extrapolation'] = True
	parameters["form_compiler"]["quadrature_degree"] = 4
	parameters["form_compiler"]["cpp_optimize"] = True
	parameters["form_compiler"]["representation"] = "uflacs"


	timer = Timer("00: forward run")
	model = problem()

	# -------UNCERTAINTYQUANTIFICATION-------------        
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

 
        # Define Joint probability density funtion
        #joint = cp.J(b_ff_dist,b_xx_dist,b_fx_dist,CC_dist,K_dist,fiber_rotation_epi_dist)
        joint = cp.J(b_ff_dist,b_xx_dist,b_fx_dist,CC_dist,K_dist,fiber_rotation_epi_dist,fiber_rotation_end_dist,sheet_rotation_epi_dist,sheet_rotation_end_dist)


        # Polynomial Chaos Expansion: three terms recursion relation
        order = 2
        Polynomials = cp.orth_ttr(order, joint)

        # Generate samples by quasi-random samples
        qmc_scheme = "S"
        nnodes = 2*len(Polynomials) #2*len(Polynomials) 
        nodes = joint.sample(nnodes, qmc_scheme)

        # Regression method in Point Collocation
        rgm_rule ="T" #Orthogonal Matching Pursuit#"LS"

	# range of pressure
        pi = 0.0
        pf = 2.
        dp = 0.05 #0.2
        pressures = np.arange(pi, pf+dp, dp)

        print "---------------------------------"
        print "PCE: parameters"
        print "Order:", order
        print "Sampling Method:", qmc_scheme
        print "Number of nodes:", nnodes
        print "Regression method:", rgm_rule
        print "Range of Pressures:", pressures
        print "---------------------------------"



        # file names
        fileroot = "UQ_PCE-PC"  + "_o" + str(order) + "_qmc" + str(qmc_scheme) + "_nnodes"+ str(nnodes) + "_rm" + str(rgm_rule) +\
        "_pi" + str(int(pi*100)) + "_pf" + str(int(pf*100)) + "_dp" + str(int(dp*100)) + \
        "_sp" + str(int(len(joint)))

        qoi_names = np.array(["dVol", "dLen", "dRen", "dRep", "dThi", "dTwi","dWV"]).T



        txt_filename = str(fileroot) + "_qoi.txt"
        time_filename = str(fileroot) + "_time.txt"
        umean_filename = str(fileroot) + "_umean.pvd"
        ustdv_filename  = str(fileroot) + "_ustd_.pvd"
        ums_filename = str(fileroot) + "_ums_.pvd"
        sffmean_filename = str(fileroot) + "_sffmean.pvd"
        sffstdv_filename  = str(fileroot) + "_sffstd_.pvd"

	sai_filename = str(fileroot) + "_SAi_.txt"
        usai1_filename =str(fileroot) + "_pressure" +  str(int(pf*100)) + "_uSAi1_.pvd"
        usait_filename =str(fileroot) + "_pressure" +  str(int(pf*100)) + "_uSAit_.pvd"
        sffsai1_filename =str(fileroot) + "_pressure" +  str(int(pf*100)) + "_sffSAi1_.pvd"
        sffsait_filename =str(fileroot) + "_pressure" +  str(int(pf*100)) + "_sffSAit_.pvd"

        hist_filename = str(fileroot) + "_pressure" +  str(int(pf*100)) + "_hist.txt"
        holdout_filename =  str(fileroot) + "_pressure" +  str(int(pf*100)) + "_hout_.txt"
        

	# open files
        txt_file = open(txt_filename,"w")
        time_file = open(time_filename, "w")
	sai_file = open(sai_filename,"w")
        hist_file = open(hist_filename,"w")
	holdout_file = open(holdout_filename,"w")

        umean_file = File(umean_filename)
        ustdv_file = File(ustdv_filename)
        ums_file = File(ums_filename)
        usai1_file = File(usai1_filename)
        usait_file = File(usait_filename)
 
        sffmean_file = File(sffmean_filename)
        sffstdv_file = File(sffstdv_filename)
        sffsai1_file = File(sffsai1_filename)
        sffsait_file = File(sffsait_filename)


        umean = Function(model.V)
        ustdv = Function(model.V)
        ums = Function(model.V)
        usai1 = Function(model.V)
        usait = Function(model.V)

        sffmean = Function(model.projection_space)
        sffstdv = Function(model.projection_space)
        sffsai1 = Function(model.projection_space)
        sffsait = Function(model.projection_space)

        # Evaluate model in sample points
        solves1 = []
        solves2 = []
        solves3 = [] 
        nnan, i= 0,0
        for s in nodes.T:
          #print "----------------"
          #print "NODE", s
          #print "----------------"
          model = problem()
          model.update_params(*s)
          pressure = pi
          while pressure <= pressures[len(pressures)-1]:
                print "----------------"
                print "PRESSURE", pressure
                print "----------------"
                qoi = model.run(pressure)
                if (i == 0)&(pressure == 0.0):

                  txt_file.write("--------------------------------- \n")
                  txt_file.write("PCE: parameters \n")
                  txt_file.write("Order: %i \n" % order)
                  txt_file.write("Sampling Method: %s \n" % qmc_scheme)
                  txt_file.write("Number of nodes: %i \n" % nnodes)
                  txt_file.write("Regression method: %s \n" % rgm_rule)
                  txt_file.write("Range of Pressures: %.2f,%.2f,%.2f \n" % (pi ,pf ,dp ))
                  txt_file.write("--------------------------------- \n")
                  txt_file.write("\n")
                  txt_file.write("Initial Values: \n" )
                  np.savetxt(txt_file, np.column_stack((qoi_names, np.array(qoi[0:7]).T.astype(np.object))) , fmt=['%s','%.8f'], delimiter =  "  ")
                  txt_file.write("\n")
                  txt_file.write("\n")

                while np.isnan(qoi[0]) == True:
                        nnan = nnan + 1
                        print "WARNING  NaN:", nnan
                        newsample = joint.sample(1, rule = "R")
                        model = problem()
                        model.update_params(*newsample.T[0])
                        pressure = pi
                        qoi = model.run(pressure)

                pressure = pressure + dp

          solves1.append(qoi[0:7])
          solves2.append(qoi[7])
          solves3.append(qoi[8])
          i = i +1


        # Create surrogate model
        #u_hat(x,t;q)=sum_n C_n(x) * P_n(q)
        qoi1_hat =  cp.fit_regression(Polynomials, nodes, solves1,rule=rgm_rule)
        qoi2_hat =  cp.fit_regression(Polynomials, nodes, solves2,rule=rgm_rule)
        qoi3_hat =  cp.fit_regression(Polynomials, nodes, solves3,rule=rgm_rule)

        #cross-validation/holdout method
        npoints = 1
        ntests = nnodes + npoints
        test_rule =  qmc_scheme
        test_samples = joint.sample(ntests, rule = test_rule)
        qoi1_samples = []
        #qoi1_hat_samples.append(qoi1_hat(*nodes))
        for ts in [test_samples.T[ntests-1]]:#[test_samples.T[ntests -1]]: #test_samples.T[ntests-npoints-1: ntests -1]:
          print ts
          model = problem()
          model.update_params(*ts)
          pressure = pi
          while pressure <= pressures[len(pressures)-1]:
                qoi = model.run(pressure)
                while np.isnan(qoi[0]) == True:
                        nnan = nnan + 1
                        print "WARNING  NaN:", nnan
                        newsample = joint.sample(1, rule = "R")
                        model = problem()
                        model.update_params(*newsample.T[0])
                        pressure = pi
                        qoi = model.run(pressure)
                        test_samples.T[ntests-1] = newsample.T[0]

                pressure = pressure + dp

          #add evaluation fmovel in point
          solves1.append(qoi[0:7])


        #metamodel evaluation in all points
        qoi1_hat_samples = []
        qoi1_hat_samples = qoi1_hat(*test_samples)
        #mse = np.array((np.array(qoi1_samples) - np.array(qoi1_hat_samples))**2).mean(axis=0).T
        #Predictive error
        pe_top = np.array((solves1 - qoi1_hat_samples.T)**2).sum(axis=0).T
        #pe_top = np.array((solves1 - qoi1_hat_samples.T)**2).sum(axis=0).T
        pe_bot = np.array((solves1 -np.mean(qoi1_hat_samples.T, axis=0))**2).sum(axis=0).T
        pe = np.divide(pe_top,pe_bot)
        #Mean Square error
        mse = (np.array((solves1[:-1] - qoi1_hat(*nodes).T)**2).sum(axis=0).T)/nnodes

        #output qoi metrics
        mean = cp.E(qoi1_hat, joint)
        var = cp.Var(qoi1_hat, joint)
        std = np.sqrt(var)
        cv = np.divide(np.array(std), np.array(mean))
        kurt = cp.Kurt(qoi1_hat,joint)
        skew = cp.Skew(qoi1_hat,joint)
        metrics = np.array([mean, std,var,cv,kurt,skew]).T
        rawdata = np.array(qoi1_hat(*joint.sample(10**5))).T#np.array(solves1)
        prodata =  np.column_stack((qoi_names, metrics.astype(np.object)))
        #prodata[:,1].astype(float)


        # saving qoi results
        txt_file.write("Pressure (kPa) = %f \n" % pressure)
        np.savetxt(txt_file, prodata, fmt=['%s','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f'], delimiter =  "  ",header = "qoi, mean, stdv, vari, cova, kurt, skew")
        txt_file.write("\n")
        txt_file.write("\n")
        #np.savetxt(txt_file, metrics,fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f'], header = "mean, stdv, vari, cova, kurt, skew")
        # saving qoi raw data
        np.savetxt(hist_file, rawdata, fmt=['%.10f','%.10f', '%.10f','%.10f','%.10f','%.10f','%.10f'], delimiter = " ", header = "dVol, dLen, dRen, dRep, dThi, dTwi, dWV")

        # saving cross-validation/holdout method 
        holdout_file.write("Polynomial order: %i \n" % order)
        holdout_file.write("Num. test samples: %i \n" % ntests)
        holdout_file.write("Test sampling rule: %s \n" % test_rule)

        holdout_file.write("\n")
        #np.savetxt(holdout_file, mse, fmt=['%.3f','%.3f', '%.3f','%.3f','%.3f','%.3f','%.3f'], delimiter = " ", header = "dVol, dLen, dRen, dRep, dThi, dTwi, dWV")
        holdout_file.write("Predictive error: \n")
        holdout_file.writelines(["%s\n" % item  for item in pe])
        txt_file.write("\n")
        holdout_file.write("Mean Square error: \n")
        holdout_file.writelines(["%s\n" % item  for item in mse])



        # output sensitivity index 
        sai1 = cp.Sens_m(qoi1_hat,joint)
        sait = cp.Sens_t(qoi1_hat,joint)
        sai2 = cp.Sens_m2(qoi1_hat,joint)

        # save sensitivity index 
        sai_file.write("Pressure (kPa) = %f \n" % pressure)
	sai_file.write("First order Sensitivity index \n" )
        #np.savetxt(sai_file, sai1.T, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f'], header = " b_ff, b_xx, b_fx, C, K, fr_epi")
        np.savetxt(sai_file, sai1.T, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f','%.5f','%.5f','%.5f'], header = " b_ff, b_xx, b_fx, C, K, fr_epi, fr_endo, sr_epi, sr_endo")

        sai_file.write("\n")
        sai_file.write("Total-effect index \n" )
	#np.savetxt(sai_file, sait.T, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f'], header = " b_ff, b_xx, b_fx, C, K, fr_epi")
        np.savetxt(sai_file, sait.T, fmt=['%.5f','%.5f', '%.5f','%.5f','%.5f','%.5f','%.5f','%.5f','%.5f'], header = " b_ff, b_xx, b_fx, C, K, fr_epi,fr_endo, sr_epi, sr_endo")

        sai_file.write("\n")
        sai_file.write("\n")

	#output displacement/stress metrics
        #sai1 = cp.Sens_m(qoi_hat[0],joint)
        u_mean = cp.E(qoi2_hat, joint)
        u_var = cp.Var(qoi2_hat, joint)
        u_stdv = np.sqrt(u_var)
	u_sai1 = cp.Sens_m(qoi2_hat,joint)
	u_sai2 = cp.Sens_m2(qoi2_hat,joint)
	u_sait = cp.Sens_t(qoi2_hat,joint)

        sff_mean = cp.E(qoi3_hat, joint)
        sff_var = cp.Var(qoi3_hat, joint)
        sff_stdv = np.sqrt(sff_var)
        sff_sai1 = cp.Sens_m(qoi3_hat,joint)
        sff_sai2 = cp.Sens_m2(qoi3_hat,joint)
        sff_sait = cp.Sens_t(qoi3_hat,joint)


        #second order index          
        qois = ["dv", "dl", "dRen", "dRep", "dT", "dTwi", "dw"]
        sai2_reshape = sai2.reshape(len(qois),int(len(joint)),int(len(joint)))

        for i in range(len(qois)):
          sai2_file = str(fileroot)  + "_SAi2_" + str(qois[i]) + ".txt"
          df = pd.DataFrame(data=sai2_reshape[i])
          df.to_csv(sai2_file, sep=' ', header=" b_ff, b_xx, b_fx, C, K, fr_epi,fr_endo, sr_epi, sr_endo",float_format='%.5f', index=True)


        # save displacement/stress results
        #umean.vector()[:] = u_mean
        #ustdv.vector()[:] = u_stdv
        #ums.vector()[:] = u_mean + u_stdv
	#umean_file << umean
        #ustdv_file << ustdv
        #ums_file << ums

        #sffmean.vector()[:] = sff_mean
        #sffstdv.vector()[:] = sff_stdv
        #sffmean_file << sffmean
        #sffstdv_file << sffstdv


	#for param in range(len(nodes)):
        #        usai1.vector()[:] = u_sai1[param]
        #        usait.vector()[:] = u_sait[param]
        #        usai1_file << usai1
        #        usait_file << usait
        #        sffsai1.vector()[:] = sff_sai1[param]
        #        sffsait.vector()[:] = sff_sait[param]
        #        sffsai1_file << sffsai1
        #        sffsait_file << sffsait






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
