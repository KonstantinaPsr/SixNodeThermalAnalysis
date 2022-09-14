

import numpy as np
import matplotlib.pyplot as plt
# GLOBAL_PARAMETERS
JD = 2459946									#calculated via https://planetcalc.com/503/ 
RI = 97.786*3.1415926/180						# Obit's inclination
Re = 6378										# Earth's radius
h = 600											# Orbit's height
R = 6978 				
J2 = 1082.62*1E-6
mu 		= 3.986004415e05    					# Gravitational constant (km3/s2)
vel = np.sqrt(mu/(R**3))#0.001081 				# Orbit's angular velocity in rad/s
v = -(3/2)*(Re**2/R**2)*J2*vel*np.cos(RI) 		#rad/s
#v = -(3/2)*(Re**2/R**2)*J2*vel*np.cos(RI)*(3600*24) 		#rad/day
omega = 30				  						# Initial RAAN in degrees
P  = 5812.0#2*3.1415926/vel   					# Orbit's Period
beta_cr = 66.06*3.1415926/180 					# Critical beta angle for specific height
kappa = 400 									#W/m2/K conservative value for Contact conductance
sigma 	= 5.67e-8								# Stefab-Boltzman constant (W/m^2 K^4)
m    	= 4000			    					#Sphere mass (g)
dt      = 1.0#1.0/vel

A_zn = 0.03 	# area in m2
A_ns = 0.03 	# area in m2
A_u = 0.01 		# area in m2
Asmall= 5e-4    # small contact area
Along = 1.5e-3  # long contact area
m1 = 1200 		# weight of big area in g
m2 = 400 		# weight of small area in g
q_sol = 1368 	#W/m2 mean value of solar radiation
qgen1 = 124.026 #W/m2 generated heat from large faces so that generated heat/ area = constant
qgen2 = 41.342  #W/m2 generated heat from small faces

#Aluminium 6061-T6
Cp    		= 0.896									# Specific heat capacity (J/g oC)
a 			= 0.96 									#absorptivity (as if coated by Aeroglaze 2306 black paint)
eps			= 0.9 									# emissivity

#Initation of Temperatures and parameters
Tzen_pr = 298.15
Tnad_pr = 298.15
Tn_pr = 298.15
Ts_pr = 298.15
Tun_pr = 298.15
Tup_pr = 298.15
Tmean = 298.15
n = 1
nn = 1
nnn = 1
nnnn = 3
m = 0
mm = 1
mmm = 3
l = 0
ll = 0
b = 1
bb = 1
y = 1
kk = 0

t = 0.0

fefile = open("Eclipse6.txt","w")
fefile.write('time 	 	eclipse		\n')
yfile = open("SixNode.txt","w")
yfile.write('time 	 	Tzen 			Tnad 			Tun 			Tup 			Tn 				Ts 				Tmean		\n')
Declfile	= open('Declination6.txt', 'w')
xfile = open("beta_angle6.txt","w")
xfile.write('time 	 	beta		\n')


#time control
counter = 0
tt = 0.0
i = 1
ii = 1

#https://en.wikipedia.org/wiki/Position_of_the_Sun
while (JD < 2460631):
	n = JD - 2451545.0 #Number of days since Greenwich noon, Terrestrial time, on 1 January 2000
#####################################
# In the Ecliptic Coordinate System #
#####################################
	L = 280.460 + 0.9856474*n # The mean longitude of the Sun, corrected for the aberration of light in degrees
	g = 357.528 + 0.9856003*n # The mean anomaly of the Sun

# Put L and g in the range 0° to 360° by adding or subtracting multiples of 360° as needed.
	while (L < 0):
		L = L+360*10
		while L>360:
			L = L-360
	while (g < 0):
		g = g+360*10
		while g>360:
			g = g-360

	L = np.deg2rad(L)
	g = np.deg2rad(g)
#the ecliptic longitude of the Sun is
	lamda = L +np.deg2rad(1.915)*np.sin(g) + np.deg2rad(0.020)*np.sin(2.0*g)
#The ecliptic latitude of the Sun is nearly:
	b = 0.0 # as the ecliptic latitude of the Sun never exceeds 0.00033°

# and the distance of the Sun from the Earth, in astronomical units, is:
	RR = np.deg2rad(1.00014) - np.deg2rad(0.01671)*np.cos(g)-np.deg2rad(0.00014)*np.cos(2.0*g)

# The obliquity of the ecliptic can be approximated:
	eps = 23.439 - 0.0000004*n
	eps = np.deg2rad(eps)
#######################################
# In the Equatorial Coordinate System # 
#######################################
	
	omega_s = np.arctan(np.cos(eps)*np.tan(lamda))
	decl_s  = np.arcsin(np.sin(eps)*np.sin(lamda))
	#decl_s = decl_s*180/3.1415926

	Declfile.write("%s %s \n" % (decl_s, omega_s))

	# RAAN = omega +0.985447*kk			#kk in days
	# RAAN = RAAN*3.1415926/180	
	# if (RAAN>2.0*3.1415926):
	# 	RAAN = RAAN - 2.0*3.1415926

#for t in range(360):

	
	while (tt<=86400*ii):
		RAAN = omega + v*tt
		RAAN = omega*3.1415926/180
		if (RAAN<0):
			RAAN = RAAN +2.0*3.1415926
	
		x = (np.cos(decl_s)*np.sin(RI)*np.sin(RAAN) - np.sin(decl_s)*np.cos(omega_s)*np.sin(RI)*np.cos(RAAN) - np.sin(decl_s)*np.sin(omega_s)*np.cos(RI))
		beta = np.arcsin(x)
		if (abs(beta)<beta_cr): 
			fe = (1/180)*np.arccos(np.sqrt(2*Re*h+h**2)/((Re+h)*np.cos(beta)))
		else:

			fe = 0.0	
		fefile.write("%s %s \n"%(tt , fe))

		if (beta<30*3.1415926/180): 
			aa = 0.14
			qIR = 228
		else:
			aa = 0.19
			qIR = 218	

	#14.8648*P):
		# Nadir Face
		if (n*P/4 < tt <nn*(P/2)*(1-fe) or nnn*(P/2)*(1+fe)<tt<nnnn*P/4):
			Fnad = -np.cos(vel*tt)*np.cos(beta)
		else:
			Fnad = 0.0
		
			
		# Zenith Face
		if (m*P<tt<mm*P/4) or (mmm*P/4< tt <(m+1)*P):
			Fzen = np.cos(vel*tt)*np.cos(beta)
			
		else:
			Fzen = 0.0
		
			

		# Positive Velocity Face 
		if (b*(P/2)*(1+fe)<tt<bb*P):
			Fup = -np.sin(vel*tt)*np.cos(beta)
		else:
			Fup = 0.0
			
		
		# Negative Velocity Face 
		if (l*P<tt<(ll+1)*(P/2)*(1-fe)):
			Fun = np.sin(vel*tt)*np.cos(beta)
		else:
			Fun = 0.0
		
			
		# North and South Face
		if (y*(P/2)*(1-fe)<tt<y*(P/2)*(1+fe)):
			Fns = np.sin(beta)
				
		else:
			Fns = 0.0
		
			
			
		if (tt==i*P):
		
			i+=1
			n+=4
			nn +=2
			nnn+=2
			nnnn +=4
			m+=1
			mm+=4
			mmm+=4
			b+=2
			bb+=1
			y +=2
			l+=1
			ll+=2
			#print('yes')	

		Qzen = Fzen*A_zn*q_sol*a - sigma*eps*A_zn*pow(Tzen_pr,4) +qgen1*A_zn  + kappa*(Along*(Ts_pr-Tzen_pr)+Along*(Tn_pr-Tzen_pr)+Asmall*(Tup_pr-Tzen_pr)+Asmall*(Tun_pr-Tzen_pr)) 
		Tzen = Tzen_pr + (dt/(Cp*m1))*Qzen
		Tzen_pr = Tzen
		
		Qnad = (aa+Fnad)*A_zn*q_sol*a +qIR*A_zn - sigma*eps*A_zn*pow(Tnad_pr,4)+qgen1*A_zn  + kappa*(Along*(Ts_pr-Tnad_pr)+Along*(Tn_pr-Tnad_pr)+Asmall*(Tup_pr-Tnad_pr)+Asmall*(Tun_pr-Tnad_pr))
		Tnad = Tnad_pr + (dt/(Cp*m1))*Qnad
		Tnad_pr = Tnad
		
		Qup = Fup*A_u*q_sol*a - sigma*eps*A_u*pow(Tup_pr,4)+qgen2*A_u  + kappa*(Along*(Tnad_pr-Tup_pr)+Along*(Tzen_pr-Tup_pr)+Along*(Tn_pr-Tup_pr)+Along*(Ts_pr-Tup_pr))
		Tup = Tup_pr + (dt/(Cp*m2))*Qup
		Tup_pr = Tup

		Qun = Fun*A_u*q_sol*a - sigma*eps*A_u*pow(Tun_pr,4)+qgen2*A_u + kappa*(Along*(Tnad_pr-Tun_pr)+Along*(Tzen_pr-Tun_pr)+Along*(Tn_pr-Tun_pr)+Along*(Ts_pr-Tun_pr))
		Tun = Tun_pr + (dt/(Cp*m2))*Qun
		Tun_pr = Tun

		Qn  = Fns*A_ns*q_sol*a - sigma*eps*A_ns*pow(Tn_pr,4)+qgen1*A_ns  + kappa*(Along*(Tnad_pr-Tn_pr)+Along*(Tzen_pr-Tn_pr)+Asmall*(Tup_pr-Tn_pr)+Asmall*(Tun_pr-Tn_pr))
		Tn = Tn_pr + (dt/(Cp*m1))*Qn
		Tn_pr = Tn

		Qs  =   - sigma*eps*A_ns*pow(Ts_pr,4)+qgen1*A_ns + kappa*(Along*(Tnad_pr-Ts_pr)+Along*(Tzen_pr-Ts_pr)+Asmall*(Tup_pr-Ts_pr)+Asmall*(Tun_pr-Ts_pr))+ Fns*A_ns*q_sol*a #+ kappa*(Along*(Tnad_pr-Ts_pr)+Along*(Tzen_pr-Ts_pr)+Asmall*(Tup_pr-Ts_pr)+Asmall*(Tun_pr-Ts_pr))
		Ts = Ts_pr + (dt/(Cp*m1))*Qs
		Ts_pr = Ts

		Tmean = (Tzen+Tnad+Tup+Tun+Tn+Ts)/6

		if counter%1000== 0 :

			yfile.write("%s 		%s 			%s 		%s 			%s 			%s 			%s 			%s\n" %(tt, Tzen, Tnad, Tun, Tup, Tn, Ts, Tmean))

			
		tt +=float(dt)
		counter += 1
	xfile.write("%s 		%s\n" %(tt, beta*180/3.1415926))

	JD+=1
	kk+=1
	ii+=1

xfile.close()
yfile.close()
Declfile.close()
fefile.close()

