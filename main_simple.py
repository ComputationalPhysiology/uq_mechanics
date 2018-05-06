#/usr/bin/python




if __name__ == "__main__":
	from dolfin import *
	import chaospy as cp
        import numpy  as np

	from ProblemSettings import *



	parameters['allow_extrapolation'] = True
	parameters["form_compiler"]["quadrature_degree"] = 4
	parameters["form_compiler"]["cpp_optimize"] = True
	parameters["form_compiler"]["representation"] = "uflacs"


	timer = Timer("00: forward run")
	model = problem()
        
        b_ff = Constant(6.6)
        b_xx = Constant(4.0)
        b_fx = Constant(2.6)
        CC   = Constant(1.1)      #kPa
        K    = Constant(10.0)     #kPa
        fiber_rotation_epi = 50.0  #kPa
        fiber_rotation_endo = 40.0
        sheet_rotation_epi = 65.0
        sheet_rotation_endo = 25.0
        alpha_noise = 0.0
        beta_noise = 0.0
        
	model.update_params(b_ff,b_xx,b_fx,CC,K, fiber_rotation_epi, fiber_rotation_endo,sheet_rotation_epi,sheet_rotation_endo,alpha_noise, beta_noise)
        
	#pressure = 0.1 #kPa
        file = File("displacement.pvd")
        displacement = Function(model.V)

        pressures = np.arange(0.0, 0.1, 0.1)
	print pressures 
	for pressure in pressures:#[0.00,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.040,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08]:
		print pressure
		qoi = model.run(pressure)
		#desplacement = Function(model.V)
        	displacement.vector()[:] = qoi[6]
#		file = File("desplacement.pvd")	 
		file << displacement
        

	time = timer.stop()
	print "time = ", time

	
	
