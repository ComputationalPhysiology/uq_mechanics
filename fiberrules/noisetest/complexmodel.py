#!/usr/bin/python




if __name__ == "__main__":
	from dolfin import *
	import chaospy as cp
        import numpy  as np

	from ProblemSettings import *

	parameters["form_compiler"]["quadrature_degree"] = 3
	parameters["form_compiler"]["cpp_optimize"] = True
	parameters["form_compiler"]["representation"] = "uflacs"


	timer = Timer("00: forward run")
	#geoparam = inputs.geometry_parameters()
	model = problem()

	# Define random distribution for each parameter
	b_ff_dist = cp.Normal(mu = 6.6, sigma = 0.99)
	b_xx_dist = cp.Normal(mu = 4.0, sigma = 0.6)
	b_fx_dist = cp.Normal(mu = 2.6, sigma = 0.39)
	CC_dist   = cp.Normal(mu = 1.1 , sigma = 0.165)#cp.Lognormal (mu = -0.049, sigma = 0.538)#(mu = 1.1 , sigma = 0.165)      #kPa
	K_dist    = cp.Normal(mu = 10.0, sigma = 1.5)#cp.Lognormal(mu = 2.2964, sigma = 0.11)#cp.Normal( mu = 10.0, sigma = 1.5)           #kPa
	fiber_rotation_epi_dist = cp.Normal(mu = 60,sigma = 9)                #deg


	b_ff = Constant(6.6)
	b_xx = Constant(4.0)
	b_fx = Constant(2.6)
	CC   = Constant(1.1)      #kPa
	K    = Constant(10.0)     #kPa
        fiber_rotation_epi = Constant(50)

        # Define Joint probability density funtion
        joint = cp.J(b_ff_dist,b_xx_dist,b_fx_dist,CC_dist,K_dist, fiber_rotation_epi_dist)

        # Monte Carlo	
        nsamples= 1
        rule = "R"
        samples = joint.sample(nsamples, rule=rule)
	#samples = [b_ff,b_xx,b_fx,CC,K,fiber_rotation_epi]
        #model.update_params(alpha,b_ff,b_xx,b_fx,CC,K)
        #u = model.run(2)
	file = File("displacement.pvd")
	qoi_mc = []
	for s in samples.T:
          #model.update_params(*s)
	  model.update_params() 
	  for pressure in [0.0,0.1]:#,0.2,0.3,0.4]:#,0.5,0.6,0.7,0.8,0.9,1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1]:
          	#file << model.run(pressure)
                qoi = model.run(pressure)
		print qoi


          #plot(u,mesh,interactive=True)
	  #qoi_mc.append(model.run(0.01))

	

        # output  metrics
        mean = np.mean(qoi_mc, 0)
        std = np.std(qoi_mc, 0)
        var = np.var(qoi_mc, 0)
        cv = np.divide(np.array(std), np.array(mean))
        error = std / np.sqrt(nsamples)


	time = timer.stop()
	print "time = ", time
