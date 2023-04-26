import numpy as np
import copy
from sympy import var, solve
from scipy.optimize import fsolve


class BaseThermodynamicModel:

  ###  Attributes
  R_universal      = 8.31446261815324      # Universal gas constant [j/(mol k)]
  R_specific       = -1.0                  # Specific gas constant [J/(kg K)]
  molecular_weight = -1.0                  # Molecular weight [kg/mol]
  gamma            = -1.0                  # Ratio of heat capacities [-]

  #####

  ### Constructor
  def __init__(self):
    pass

  ### Methods
  def calculatePressureFromTemperatureDensity(self, T, rho):
    
    return 0.0


  def calculateTemperatureFromPressureDensity(self, P, rho):

    return 0.0
  

  def calculateTemperatureFromPressureDensityWithInitialGuess(self, T, P, rho):

    return 0.0

  
  def calculateInternalEnergyFromPressureTemperatureDensity(self, P, T, rho):

    return 0.0


  def calculatePressureTemperatureFromDensityInternalEnergy(self, P, T, rho, e):

    return 0.0


  def calculateDensityInternalEnergyFromPressureTemperature(self, rho, e, P, T):

    return 0.0


  def calculateSpecificHeatCapacities(self, c_v, c_p, P, T, rho):

    return 0.0 


  def calculateHeatCapacitiesRatio(self, P, T, rho):

    return 0.0 


  def calculateSoundSpeed(self, P, T, rho):

    return 0.0


  def calculateVolumeExpansivity(self, T, bar_v):

    return 0.0


  def calculateIsothermalCompressibility(self, T, bar_v):

    return 0.0

  
  def calculateIsentropicCompressibility(self, P, T, bar_v):
  
    return 0.0

###################################################################
##################### IDEAL gas ###################################
###################################################################

class IdealGasModel(BaseThermodynamicModel):

  ###  Attributes

  ### Constructor
  def __init__(self, R_specific, gamma):
    super(BaseThermodynamicModel,self).__init__()  
    self.R_specific       = R_specific
    self.gamma            = gamma    
    self.molecular_weight = self.R_universal/self.R_specific

  ### Methods
  def calculatePressureFromTemperatureDensity(self, T, rho):

    # Equation of state
    P = rho*self.R_specific*T

    return P


  def calculateTemperatureFromPressureDensity(self, P, rho):
    
    # Equation of state
    T = P/(self.R_specific*rho)

    return T
  

  def calculateTemperatureFromPressureDensityWithInitialGuess(self, T, P, rho):

    # Equation of state
    T = P/(self.R_specific*rho)

    return T


  def calculateInternalEnergyFromPressureTemperatureDensity(self, P, T, rho):

    #Equation of Specific heat at constant volume
    c_v = self.R_specific/(self.gamma - 1.0)

    #Equation of Internal Energy 
    e = c_v*T

    return e


  def calculatePressureTemperatureFromDensityInternalEnergy(self, P, T, rho, e):

    c_v = self.R_specific/(self.gamma - 1.0)
    P = e*rho*(self.gamma - 1.0 )
    T = e/c_v

    return P, T


  def calculateDensityInternalEnergyFromPressureTemperature( self, rho, e, P, T):

    c_v = self.R_specific/(self.gamma - 1.0)

    e   = c_v*T
    rho = P/(e*(self.gamma - 1.0))

    return rho, e 


  def calculateSpecificHeatCapacities(self, c_v, c_p, P, T, rho):
  
    c_v = self.R_specific/(self.gamma - 1.0)
    c_p = c_v*self.gamma

    return c_v, c_p


  def calculateHeatCapacitiesRatio(self, P, rho):

    return self.gamma


  def calculateSoundSpeed(self, P, T, rho):

    sos = np.sqrt(self.gamma*P/(rho))

    return sos


  def calculateVolumeExpansivity(self, T, bar_v):

    expansivity = 1.0/T

    return expansivity


  def calculateIsothermalCompressibility(self, T, bar_v):

    isothermal_compressibility = bar_v/(self.molecular_weight*self.R_specific*T)

    return isothermal_compressibility


  def calculateIsentropicCompressibility(self, P, T, bar_v):

    dP_dT_const_v = self.molecular_weight*self.R_specific/bar_v
    dP_dv_const_T = (-1.0)*self.molecular_weight*self.R_specific*T/(bar_v*bar_v)

    isothermal_compressibility = ( -1.0 )/(bar_v*dP_dv_const_T)
    expansivity                = ( -1.0 )*(dP_dT_const_v/(bar_v*dP_dv_const_T))

    c_v = self.R_specific/(self.gamma - 1.0 )
    c_p = c_v*self.gamma

    bar_c_p = self.molecular_weight*c_p

    isentropic_compressibility = (isothermal_compressibility - ((bar_v*T*expansivity**2.0)/bar_c_p))

    return isentropic_compressibility
    

###################################################################
##################### REAL gas ####################################
###################################################################

class PengRobinsonModel(BaseThermodynamicModel):

  ### Atributes (Variables that don't change with temperature)
  acentric_factor               = -1.0   
  critical_temperature          = -1.0
  critical_pressure             = -1.0  
  critical_molar_volume         = -1.0
  NASA_coefficients             = (-1.0)*np.ones(15)
  eos_b                         = -1.0 
  eos_ac                        = -1.0
  eos_kappa                     = -1.0
  max_aitken_iter               = 1000
  aitken_relative_tolerance     = 1.0e-5
  ### Nonlinear solver parameters
  xtol   = 1.0e-10		# Relative error between two consecutive iterates
  epsfcn = 1.0e-5			# Step length for the forward-difference approximation of the Jacobian
  factor = 1.0e-1			# Parameter determining the initial step bound. Should be in the interval [0.1:100]

  ### Constructor

  def __init__(self, W, acentric_factor, critical_temperature, critical_pressure, critical_molar_volume, NASA_coefficients):
    super(BaseThermodynamicModel,self).__init__()  
    self.molecular_weight           = W
    self.R_specific                 = self.R_universal/W     
    self.acentric_factor            = acentric_factor
    self.critical_temperature       = critical_temperature
    self.critical_pressure          = critical_pressure
    self.critical_molar_volume      = critical_molar_volume
    self.NASA_coefficients          = NASA_coefficients

    self.eos_b  = 0.077796*( self.R_universal*self.critical_temperature/self.critical_pressure )
    #print('eos_b: {}'.format(eos_b))

    self.eos_ac = 0.457236*(( self.R_universal*self.critical_temperature)**2.0/(self.critical_pressure))
    #print('eos_ac: {}'.format(eos_ac))

    if acentric_factor > 0.49:
      self.eos_kappa = 0.379642 + 1.48503*self.acentric_factor - 0.164423*self.acentric_factor**2.0 + 0.016666*self.acentric_factor**3.0 
    else:
      self.eos_kappa = 0.37464 + 1.54226*self.acentric_factor - 0.26992*self.acentric_factor**2.0


  ### Methods

  def calculateTemperatureFromPressureDensity( self, P, rho ):

    # Calculate molar volume 
    bar_v = self.molecular_weight/rho
    
    # Calculate temperature guess using ideal-gas model
    T = P*bar_v/self.R_universal

    x_0 = T

    for iter in range(self.max_aitken_iter):
      x_1 = ((bar_v - self.eos_b) / self.R_universal) * (P + (self.calculate_eos_a(x_0) /((bar_v**2.0) + 2.0 * self.eos_b * bar_v - (self.eos_b**2.0))))
      x_2 = ((bar_v - self.eos_b) / self.R_universal) * (P + (self.calculate_eos_a(x_1) /((bar_v**2.0) + 2.0 * self.eos_b * bar_v - (self.eos_b**2.0))))
    
      denominator = x_2 - 2.0 * x_1 + x_0

      #T = x_2 - (pow(x_2 - x_1, 2.0)) / denominator
      T = x_2 - (x_2 - x_1)**2.0 / (denominator + 1.0e-10)
    
      if abs((T - x_2) / T) < self.aitken_relative_tolerance:
          break  # If the result is within tolerance, leave the loop!
          
      x_0 = T  # Otherwise, update x_0 to iterate again...
      
    return T


  def calculateTemperatureFromPressureDensityWithInitialGuess( self, P, rho ):

    # Calculate molar volume 
    bar_v = self.molecular_weight/rho
    
    # Calculate temperature guess using ideal-gas model
    T = P*self.bar_v/self.R_universal

    x_0 = T

    for iter in range(self.max_aitken_iter):
      x_1 = ((bar_v - self.eos_b) / self.R_universal) * (P + (self.calculate_eos_a(x_0) /((bar_v**2.0) + 2.0 * self.eos_b * bar_v - (self.eos_b**2.0))))
      x_2 = ((bar_v - self.eos_b) / self.R_universal) * (P + (self.calculate_eos_a(x_1) /((bar_v**2.0) + 2.0 * self.eos_b * bar_v - (self.eos_b**2.0))))
    
      denominator = x_2 - 2.0 * x_1 + x_0

      #T = x_2 - (pow(x_2 - x_1, 2.0)) / denominator
      T = x_2 - (x_2 - x_1)**2.0 / (denominator + 1.0e-10)
    
      if abs((T - x_2) / T) < self.aitken_relative_tolerance:
          break  # If the result is within tolerance, leave the loop!

      x_0 = T  # Otherwise, update x_0 to iterate again...
    
    return T


  def calculateInternalEnergyFromPressureTemperatureDensity(self, P, T, rho):

    bar_v = self.molecular_weight/rho
    e = (1.0/self.molecular_weight)*self.calculateMolarInternalEnergyFromPressureTemperatureMolarVolume(P, T, bar_v)

    return e

  def n_solve( self, functions, variables, norm_factors ):
    func = lambda x : [ f(*x) for f in functions ]
    variables = fsolve( func, variables, xtol = self.xtol, epsfcn = self.epsfcn, factor = self.factor, diag = norm_factors )

    return variables

  def calculatePressureTemperatureFromDensityInternalEnergy(self, P, T, rho, e):

    bar_v = self.molecular_weight/rho

    ### Define functions
    functions = []
    functions.append( lambda variable_P, variable_T : ( ( self.calculatePressureFromTemperatureDensity(variable_T,rho) ) - variable_P )/variable_P )
    functions.append( lambda variable_P, variable_T : ( ( self.calculateMolarInternalEnergyFromPressureTemperatureMolarVolume(variable_P,variable_T,bar_v)/self.molecular_weight ) - e )/e )

    ### Initialize variables: P & T
    variables = np.zeros( 2 )
    variables[0] = P    # Use input P value as initial guess
    variables[1] = T    # Use input T value as initial guess

    ### Set normalization factors of Jacobian's diagonal: P & T
    norm_factors = np.zeros( 2 )
    norm_factors[0] = copy.deepcopy( abs( variables[0] ) )
    norm_factors[1] = copy.deepcopy( abs( variables[1] ) )

    ### Solve nonlinear system
    variables = self.n_solve( functions, variables, norm_factors )
    #print( variables )

    ### Assign solution to P & T
    P = variables[0]
    T = variables[1]
    #print( P, T)

    return P, T


  def calculateDensityInternalEnergyFromPressureTemperature(self, P, T):

    # Auxiliar parameters
    eos_a  = self.calculate_eos_a(T)
    eos_en = P*self.eos_b + self.R_universal*T
    a      = P
    b      = 2.0*P*self.eos_b - eos_en
    c      = (-1.0)*P*self.eos_b*self.eos_b - 2.0*eos_en*self.eos_b + eos_a
    d      = (-1.0)*self.eos_b*(eos_a + (-1.0)*eos_en*self.eos_b)

    print(f"a = {a}\n")
    print(f"b = {b}\n")
    print(f"c = {c}\n")
    print(f"d = {d}\n")
    
    print(f"self.eos_b = {self.eos_b}\n")
    print(f"eos_a = {eos_a}\n")
    print(f"eos_en = {eos_en}\n")

    # Cubic solve to calculate rho
    v_1, v_2, v_3 = self.calculateRootsCubicPolynomial(a, b, c, d)
    rho = self.molecular_weight/(v_1.real)     # First root is always real

    # Calculate e
    bar_v = self.molecular_weight/rho
    e = (1.0/self.molecular_weight)*self.calculateMolarInternalEnergyFromPressureTemperatureMolarVolume(P, T, bar_v)
    
    return rho, e


  def calculateSpecificHeatCapacities(self, P, T, rho):

    bar_v = self.molecular_weight/rho
    std_bar_c_p = self.calculateMolarStdCpFromNASApolynomials(T)
    std_bar_c_v = std_bar_c_p - self.R_universal
    
    c_v = (1.0/self.molecular_weight)*(std_bar_c_v + self.calculateDepartureFunctionMolarCv(P, T, bar_v))
    c_p = (1.0/self.molecular_weight)*(std_bar_c_p + self.calculateDepartureFunctionMolarCp(P, T, bar_v))

    return c_v, c_p


  def calculateHeatCapacitiesRatio(self, P, rho):
    
    bar_v = self.molecular_weight/rho
    T = self.calculateTemperatureFromPressureDensity(P, rho)
    std_bar_c_p = self.calculateMolarStdCpFromNASApolynomials(T)
    std_bar_c_v = std_bar_c_p - self.R_universal

    c_v = (1.0/self.molecular_weight) * (std_bar_c_v + self.calculateDepartureFunctionMolarCv(P, T, bar_v))
    c_p = (1.0/self.molecular_weight) * (std_bar_c_p + self.calculateDepartureFunctionMolarCp(P, T, bar_v))

    gamma = c_p/c_v
    
    return gamma


  def calculateSoundSpeed(self, P, T, rho):
    
    bar_v = self.molecular_weight/rho

    sos = np.sqrt( 1.0/(rho*self.calculateIsentropicCompressibility(P, T, bar_v)) )

    return sos
  

  def calculateVolumeExpansivity(self, T,bar_v ):

    dP_dT_const_v = self.calculateDPDTConstantMolarVolume(T,bar_v)
    dP_dv_const_T  = self.calculateDPDvConstantTemperature(T,bar_v)

    Expansivity = (-1.0)*(dP_dT_const_v/(bar_v*dP_dv_const_T))

    return Expansivity
  

  def calculateIsothermalCompressibility(self, T, bar_v):

    dP_dv_const_T = self.calculateDPDvConstantTemperature(T, bar_v)
    isothermal_compressibility = (-1.0)/(bar_v*dP_dv_const_T)

    return isothermal_compressibility
  

  def calculateIsentropicCompressibility(self, P, T, bar_v):

    isothermal_compressibility = self.calculateIsothermalCompressibility(T, bar_v)
    expansivity                = self.calculateVolumeExpansivity(T, bar_v)
    bar_c_p                    = self.calculateMolarStdCpFromNASApolynomials(T) + self.calculateDepartureFunctionMolarCp(P, T, bar_v)
      
    isentropic_compressibility = (isothermal_compressibility - ((bar_v * T * (expansivity ** 2.0)) / bar_c_p))
    
    return isentropic_compressibility
  

  def calculatePressureFromTemperatureDensity(self, T, rho):

    bar_v = self.molecular_weight/rho
    P = (self.R_universal*T/(bar_v - self.eos_b)) - (self.calculate_eos_a( T )/(bar_v*bar_v + 2.0*self.eos_b*bar_v - self.eos_b*self.eos_b))
    
    return P
  

  def calculateMolarInternalEnergyFromPressureTemperatureMolarVolume(self,P, T, bar_v):

    bar_e = self.calculateMolarStdEnthalpyFromNASApolynomials(T) + self.calculateDepartureFunctionMolarEnthalpy(P, T, bar_v) - P * bar_v
    
    return bar_e
    

  def calculate_eos_a( self, T ):
     
     eos_a = (0.457*((self.R_universal*self.critical_temperature)**2)/(self.critical_pressure))*(1+self.eos_kappa*(1-np.sqrt(T/self.critical_temperature)))**2
     
     return eos_a


  def calculate_eos_a_first_derivative(self, T):

    eos_a_first_derivative = self.eos_kappa*self.eos_ac*( ( self.eos_kappa/self.critical_temperature ) - ( ( 1.0 + self.eos_kappa )/np.sqrt(T*self.critical_temperature ) ) )
    # eos_a_first_derivative = eos_kappa*eos_ac*( ( eos_kappa/critical_temperature ) - ( ( 1.0 + eos_kappa )/sqrt( T*critical_temperature ) ) );
   
    return eos_a_first_derivative


  def calculate_eos_a_second_derivative(self, T):
 
    eos_a_second_derivative = (self.eos_kappa*self.eos_ac*(1.0 + self.eos_kappa))/(2.0*np.sqrt(T**3.0)*self.critical_temperature)

    return eos_a_second_derivative


  def calculate_Z(self, P, T,bar_v):

    Z = (P*bar_v)/(self.R_universal*T)

    return Z
  

  def calculate_A(self, P, T ):

    eos_a = self.calculate_eos_a( T )

    A = (eos_a*P)/((self.R_universal*T)**2.0)

    return A
  

  def calculate_B(self, P, T):
    
    B = (self.eos_b*P)/(self.R_universal*T)

    return B
  

  def calculate_M(self, Z, B):

    M = (Z**2+2.0*B*Z-B**2)/(Z - B)
    
    return M
  

  def calculate_N(self, eos_a_first_derivative, B):

    N = eos_a_first_derivative * (B/(self.eos_b*self.R_universal))

    return N
  

  def calculateMolarStdCpFromNASApolynomials(self, T):                          #### Acabar de definir els parametres del Nasa_coefficients 
    std_bar_c_p = 0.0
    if 200.0 <= T < 1000.0:
        std_bar_c_p = self.R_universal*(self.NASA_coefficients[7] + self.NASA_coefficients[8]*T + self.NASA_coefficients[9]*T**2.0 + self.NASA_coefficients[10]*T**3.0 + self.NASA_coefficients[11]* T**4.0)
    elif 1000.0 <= T < 6000.0:
        std_bar_c_p = self.R_universal*(self.NASA_coefficients[0] + self.NASA_coefficients[1]*T + self.NASA_coefficients[2]*T**2.0 + self.NASA_coefficients[3]*T**3.0 + self.NASA_coefficients[4]* T**4.0)
    elif T < 200:
        # Assume constant temperature below T = 200 K	    
        T_min = 200.0	    
        std_bar_c_p = self.R_universal*(self.NASA_coefficients[7] + self.NASA_coefficients[8]*T_min + self.NASA_coefficients[9]*T_min**2.0 + self.NASA_coefficients[10]*T_min**3.0 + self.NASA_coefficients[11]*T_min**4.0)
    else:
        print(f"\nNASA 7-coefficient polynomials for std bar c_p. T = {T} is above 6000 K.\n\n")
        exit()
    return std_bar_c_p
      
  def calculateMolarStdEnthalpyFromNASApolynomials(self, T):
    std_bar_h = 0.0

    if T >= 200.0 and T < 1000.0:
      #std_bar_h = R_universal*T*( NASA_coefficients[7] + NASA_coefficients[8]*T/2.0 + NASA_coefficients[9]*pow( T, 2.0 )/3.0 + NASA_coefficients[10]*pow( T, 3.0 )/4.0 + NASA_coefficients[11]*pow( T, 4.0 )/5.0 + NASA_coefficients[12]/T) - R_universal*NASA_coefficients[14];
      std_bar_h = self.R_universal*T*( self.NASA_coefficients[7] + self.NASA_coefficients[8]*T/2.0 + self.NASA_coefficients[9]*(T**2.0)/3.0 + self.NASA_coefficients[10]*(T**3.0)/4.0 + self.NASA_coefficients[11]*(T**4.0)/5.0 + self.NASA_coefficients[12]/T )
    elif T >= 1000.0 and T < 6000.0:
      #//std_bar_h = R_universal*T*( NASA_coefficients[0] + NASA_coefficients[1]*T/2.0 + NASA_coefficients[2]*pow( T, 2.0 )/3.0 + NASA_coefficients[3]*pow( T, 3.0 )/4.0 + NASA_coefficients[4]*pow( T, 4.0 )/5.0 + NASA_coefficients[5]/T ) - R_universal*NASA_coefficients[14];
      std_bar_h = self.R_universal*T*( self.NASA_coefficients[0] + self.NASA_coefficients[1]*T/2.0 + self.NASA_coefficients[2]*(T**2.0)/3.0 + self.NASA_coefficients[3]*(T**3.0)/4.0 + self.NASA_coefficients[4]*(T**4.0)/5.0 + self.NASA_coefficients[5]/T )
    elif T < 200.0:
      T_min = 200.0

      #std_bar_h_min   = R_universal*T_min*( NASA_coefficients[7] + NASA_coefficients[8]*T_min/2.0 + NASA_coefficients[9]*pow( T_min, 2.0 )/3.0 + NASA_coefficients[10]*pow( T_min, 3.0 )/4.0 + NASA_coefficients[11]*pow( T_min, 4.0 )/5.0 + NASA_coefficients[12]/T_min ) - R_universal*NASA_coefficients[14];
      std_bar_h_min   = self.R_universal*T_min*( self.NASA_coefficients[7] + self.NASA_coefficients[8]*T_min/2.0 + self.NASA_coefficients[9]*(T_min**2.0)/3.0 + self.NASA_coefficients[10]*(T_min**3.0)/4.0 + self.NASA_coefficients[11]*(T_min**4.0)/5.0 + self.NASA_coefficients[12]/T_min )
      std_bar_h_slope = self.R_universal*( self.NASA_coefficients[7] + self.NASA_coefficients[8]*T_min + self.NASA_coefficients[9]*(T_min**2.0) + self.NASA_coefficients[10]*(T_min**3.0) + self.NASA_coefficients[11]*T_min**4.0)
      std_bar_h = std_bar_h_min + std_bar_h_slope*(T-T_min)
      
    else:
       
       print(f"\nNASA 7-coefficient polynomials for std bar c_p. T = {T} is above 6000 K.\n\n")
       exit()

    return std_bar_h


  def calculateDepartureFunctionMolarCp(self, P, T, bar_v):
    # Peng-Robinson model:
    # D.Y. Peng, D. B. Robinson 
    # A new two-constants equation of state
    # Industrial and Engineering Chemistry: Fundamental , 15 , 59-64 , 1976.

    eos_a_first_derivative  = self.calculate_eos_a_first_derivative( T )
    eos_a_second_derivative = self.calculate_eos_a_second_derivative ( T )
    Z                       = self.calculate_Z(P, T, bar_v)
    A                       = self.calculate_A(P, T)
    B                       = self.calculate_B(P, T)
    M                       = self.calculate_M(Z, B)
    N                       = self.calculate_N(eos_a_first_derivative, B)

    Delta_bar_c_p = ((self.R_universal*(M - N)**2))/((M**2.0) - 2.0*A* (Z + B)) - ((T*eos_a_second_derivative)/(2.0*np.sqrt(2.0)*self.eos_b))*np.log((Z + (1.0 - np.sqrt(2.0))*B)/(Z + (1.0 + np.sqrt(2.0))*B)) - self.R_universal
    
    return Delta_bar_c_p


  def calculateDepartureFunctionMolarCv(self, P, T, bar_v):

    eos_a_second_derivative = self.calculate_eos_a_second_derivative( T )
    Z                       = self.calculate_Z(P, T, bar_v)
    B                       = self.calculate_B(P, T)

    Delta_bar_c_v = (-1.0)*((T*eos_a_second_derivative)/(2.0*np.sqrt( 2.0)*self.eos_b))*np.log((Z+(1.0 - np.sqrt( 2.0 ))*B)/(Z+(1.0 + np.sqrt(2.0))*B))
    
    return Delta_bar_c_v


  def calculateDepartureFunctionMolarEnthalpy(self, P, T, bar_v):
    # Peng-Robinson model 
    # D. Y. Peng, D. B. Robinson 
    # A new two-constant equations of State
    # Industrial and engineering Chemistry: Fundamental, 15, 59-64 , 1976.
    eos_a                    = self.calculate_eos_a( T )
    eos_a_first_derivative   = self.calculate_eos_a_first_derivative( T )
    Z                        = self.calculate_Z(P, T, bar_v)
    B                        = self.calculate_B(P, T)

    Delta_bar_h = self.R_universal*T*(Z - 1.0) + (( eos_a - eos_a_first_derivative*T)/(2.0*np.sqrt(2.0)*self.eos_b))*np.log((Z + (1.0 - np.sqrt(2.0))*B)/(Z + (1.0 + np.sqrt(2.0))*B))
    
    return Delta_bar_h
  

  def calculateTemperatureFromPressureMolarVolume(self, P, bar_v):

    # Initial temperature guess using ideal-gas model
    T = P*bar_v/self.R_universal
        
    # Aitken’s delta-squared process:
    x_0 = T

    for iter in range(self.max_aitken_iter):
      x_1 = ((bar_v - self.eos_b) / self.R_universal) * (P + (self.calculate_eos_a(x_0) / (bar_v**2 + 2*self.eos_b*bar_v - self.eos_b**2)))
      x_2 = ((bar_v - self.eos_b) / self.R_universal) * (P + (self.calculate_eos_a(x_1) / (bar_v**2 + 2*self.eos_b*bar_v - self.eos_b**2)))

      denominator = x_2 - 2*x_1 + x_0
      T = x_2 - ((x_2 - x_1)**2) / denominator

      if abs((T - x_2) / T) < self.aitken_relative_tolerance:
        break	# If the result is within tolerance, leave the loop!
            
        x_0 = T	# Otherwise, update x_0 to iterate again ...

    return T

  
  def calculateDPDTConstantMolarVolume(self, T, bar_v):

    eos_a_first_derivative  = self.calculate_eos_a_first_derivative( T )

    dP_dT_const_v = (self.R_universal/(bar_v - self.eos_b))-(eos_a_first_derivative/(bar_v*bar_v + 2.0*bar_v*self.eos_b - self.eos_b*self.eos_b))

    return dP_dT_const_v


  def calculateDPDvConstantTemperature(self, T, bar_v):

    eos_a = self.calculate_eos_a( T )

    dP_dv_const_T = (-1.0)*((self.R_universal*T)/(bar_v-self.eos_b)**2) + (eos_a*(2.0*bar_v + 2.0*self.eos_b))/((bar_v**2.0) + 2.0*bar_v*self.eos_b - self.eos_b**2.0)**2.0
    
    return dP_dv_const_T


  def calculateRootsCubicPolynomial( self, a, b, c, d):
    
    if a == 0:
      print("The coefficient of the cube of x is 0. Please use the utility for a SECOND degree quadratic. No further action taken.")
      return None, None, None                                                   # To define a null number 

    if d == 0:
      print("One root is 0. Now divide through by x and use the utility for a SECOND degree quadratic to solve the resulting equation for the other two roots. No further action taken.")
      return None, None, None                                                   # To define a null number

    b /= a                                                                      # normalize all coeficients value of the cubic equation                                                                            
    c /= a                                                                      # normalize all coeficients value of the cubic equation 
    d /= a                                                                      # normalize all coeficients value of the cubic equation 

    disc, q, r, dum1, s, t, term1, r13 = 0, 0, 0, 0, 0, 0, 0, 0                 #intermediate parameters of the roots
    q = (3.0*c - (b*b))/9.0
    r = -(27.0*d) + b*(9.0*c - 2.0*(b*b))
    r /= 54.0
    disc = q*q*q + r*r

    term1 = b/3.0

    if disc > 0:
      s = r + np.sqrt(disc)
      if s < 0:
        s = (-1.0)*( (-s)**(1.0/3.0) )
      else:
        s = s**(1.0/3.0)
      t = r - np.sqrt(disc)
      if t < 0:
        t = (-1.0)*( (-t)**(1.0/3.0) ) 
      else:
        t = t**(1.0/3.0) 
      root_1 = complex( (-1.0)*term1 + s + t, 0.0) 
      term1 += (s + t)/2.0
      term1  = np.sqrt(3.0)*(-t+s)/2
      root_2 = complex( (-1.0)*( b/3.0 + (s + t)/2.0 ), term1 )
      root_3 = complex( (-1.0)*( b/3.0 + (s + t)/2.0 ), (-1.0)*term1 )
      #print(root_1.real)
      #print(root_2.real)
      #print(root_3.real)
      return root_1, root_2, root_3
      
    # End if (disc > 0)
    print()
    print("The PengRobinsonModel::calculateRootsCubicPolynomial has found more than one real root.")
    print("Vapor-liquid equilibrium conditions must be solved.")
    print("Another option is to avoid calculating rho from P and T.")
    print()

    root_1 = complex( 0.0, 0.0 )
    root_2 = complex( 0.0, 0.0 )
    root_3 = complex( 0.0, 0.0 )

    return root_1, root_2, root_3


################################################################################
######################### BASE TRANSPORT COEFICIENTS ###########################
################################################################################


class BaseTransportCoefficients:				#### Base transport coefficients

  ### Atributes (Variables that don't change with temperature)
  R_universal      = 8.31446261815324      # Universal gas constant [j/(mol k)]
  #R_specific       = -1.0                  # Specific gas constant [J/(kg K)]
  #molecular_weight = -1.0                  # Molecular weight [kg/mol]
  mu_value         = -1.0                  # Dynamic viscosity [Pa·s]
  kappa_value      = -1.0                  # Thermal Conductivity [ W/(m·K)]

  ### Constructor 

  def __init__(self, R_universal):
    self.R_universal           = R_universal 


  ### Methods 

  def calculateDynamicViscosity(self, P, T, rho):

    return 0.0

  def calculateThermalConductivity(self, P, T, rho): 

    return 0.0 


class ConstantTransportCoefficients(BaseTransportCoefficients):                 ### Constant transport coefficients

  ### Constructor
  def __init__(self, mu, kappa):
    super(BaseTransportCoefficients,self).__init__()
    self.mu       = mu 
    self.kappa    = kappa 

  

  ### Methods

  def calculateDynamicViscosity(self, P, T, rho):
    
    return( self.mu )

  def calculateThermalConductivity(self, P, T, rho):
    
    return( self.kappa )
  

class LowPressureGasTransportCoefficients(BaseTransportCoefficients):           ### Low-pressure gas variable transport coefficients

  ### Atributes (Variables that don't change with temperature)
  mu_0        = -1.0
  kappa_0     = -1.0 
  T_0         = -1.0 
  S_mu        = -1.0
  S_kappa     = -1.0 

  ### Constructor 

  def __init__(self, mu_0, kappa_0, T_0, S_mu, S_kappa):
    super(BaseTransportCoefficients,self).__init__()  
    self.mu_0    =  mu_0
    self.kappa_0 = kappa_0
    self.T_0     = T_0
    self.S_mu    = S_mu 
    self.S_kappa = S_kappa

  def calculateDynamicViscosity(self, P, T, rho):

    return ((self.mu_0*(T/self.T_0)**1.5)*((self.T_0 + self.S_mu)/(T + self.S_mu)))

  
  def calculateThermalConductivity(self, P, T, rho):

    return (self.kappa_0*((T/self.T_0)**1.5))*((self.T_0 + self.S_kappa)/(T + self.S_kappa))


class HighPressureTransportCoeficients(BaseTransportCoefficients):              ### High-pressure transport coefficients

  ### Atributes (Variables that don't change with temperature)

  molecular_weight                = -1.0 
  critical_temperature            = -1.0
  critical_molar_volume           = -1.0
  acentric_factor                 = -1.0
  dipole_moment                   = -1.0
  association_factor              = -1.0
  NASA_coefficients               = (-1.0)*np.ones(15)

  ### Constructor

  def __init__(self, molecular_weight, acentric_factor, critical_temperature, critical_molar_volume, NASA_coefficients, dipole_moment,association_factor):
    super(BaseTransportCoefficients,self).__init__()  
    self.molecular_weight           = molecular_weight
    self.critical_temperature       = critical_temperature
    self.critical_molar_volume      = critical_molar_volume   
    self.acentric_factor            = acentric_factor
    self.dipole_moment              = dipole_moment
    self.association_factor         = association_factor  
    self.NASA_coefficients          = NASA_coefficients
    
    
    ### Adimensional dipole moment -- Poling et al. The properties of gases and liquids. McGraw-Hill, 2001.
    self.adimensional_dipole_moment = 131.3*( self.dipole_moment/np.sqrt( ( 1.0e6*self.critical_molar_volume )*self.critical_temperature ) )
      
    ### Viscosity mu -- Poling et al. The properties of gases and liquids. McGraw-Hill, 2001. (9.40, Table 9-6)
    a1_mu  = 6.324;   a2_mu  = 1.210e-3;    a3_mu  = 5.283;   a4_mu  = 6.623;   a5_mu  = 19.745    
    a6_mu  = -1.900;  a7_mu  = 24.275;      a8_mu  = 0.7972;  a9_mu  = -0.2382; a10_mu = 0.06863   
    ####### Intenta posar-ho així perque 
    b1_mu  = 50.412;  b2_mu  = -1.154e-3 ; b3_mu  = 254.209;  b4_mu  = 38.096 ;  b5_mu  = 7.630;   b6_mu  = -12.537; b7_mu  = 3.450;  b8_mu  = 1.117;   b9_mu  = 0.06770; b10_mu = 0.3479; c1_mu  = -51.680; c2_mu  = -6.257e-3;  c3_mu  = -168.48; c4_mu  = -8.464; d10_mu = -0.727;
    c5_mu  = -14.354; c6_mu  = 4.985     ; c7_mu  = -11.291;  c8_mu  = 0.01235;  c9_mu  = -0.8163; c10_mu = 0.5926 ; d1_mu  = 1189.0; d2_mu  = 0.03728; d3_mu  = 3898.0;  d4_mu = 31.42;   d5_mu  = 31.53;   d6_mu  = -18.15;     d7_mu  = 69.35;   d8_mu  = -4.117; d9_mu  = 4.025;

    self.E1_mu  = a1_mu  + b1_mu*self.acentric_factor  + c1_mu*self.adimensional_dipole_moment**4.0  + d1_mu*self.association_factor
    self.E2_mu  = a2_mu  + b2_mu*self.acentric_factor  + c2_mu*self.adimensional_dipole_moment**4.0  + d2_mu*self.association_factor
    self.E3_mu  = a3_mu  + b3_mu*self.acentric_factor  + c3_mu*self.adimensional_dipole_moment**4.0  + d3_mu*self.association_factor
    self.E4_mu  = a4_mu  + b4_mu*self.acentric_factor  + c4_mu*self.adimensional_dipole_moment**4.0  + d4_mu*self.association_factor
    self.E5_mu  = a5_mu  + b5_mu*self.acentric_factor  + c5_mu*self.adimensional_dipole_moment**4.0  + d5_mu*self.association_factor
    self.E6_mu  = a6_mu  + b6_mu*self.acentric_factor  + c6_mu*self.adimensional_dipole_moment**4.0  + d6_mu*self.association_factor
    self.E7_mu  = a7_mu  + b7_mu*self.acentric_factor  + c7_mu*self.adimensional_dipole_moment**4.0  + d7_mu*self.association_factor
    self.E8_mu  = a8_mu  + b8_mu*self.acentric_factor  + c8_mu*self.adimensional_dipole_moment**4.0  + d8_mu*self.association_factor
    self.E9_mu  = a9_mu  + b9_mu*self.acentric_factor  + c9_mu*self.adimensional_dipole_moment**4.0  + d9_mu*self.association_factor
    self.E10_mu = a10_mu + b10_mu*self.acentric_factor + c10_mu*self.adimensional_dipole_moment**4.0 + d10_mu*self.association_factor

    ### Thermal conductivity k -- Poling et al. The properties of gases and liquids. McGraw-Hill, 2001. (10.23, Table 10-3)
    a1_k = 2.4166*1.00;  a3_k = 6.6107*1.00; a5_k = 7.9274*0.10; a7_k = 9.1089*10.0;  b2_k = -1.5094*1.00; b4_k = -8.9139*1.00; b6_k = 1.2801*10.00; c1_k = -9.1858*0.10; c3_k = 6.4760*10.00; c5_k = -6.9369*0.10; c7_k = -5.4217*10.0; d2_k = 6.9983*10.00; d4_k = 7.4344*10.00; d6_k = 6.5529*10.00
    a2_k = -5.0924*0.1;  a4_k = 1.4543*10.0; a6_k = -5.8634*1.0; b1_k = 7.4824*0.100; b3_k = 5.6207*1.000; b5_k = 8.2019*0.100; b7_k = 1.2811*100.0; c2_k = -4.9991*10.0; c4_k = -5.6379*1.00; c6_k = 9.5893*1.000; d1_k = 1.2172*100.0; d3_k = 2.7039*10.00; d5_k = 6.3173*1.000; d7_k = 5.2381*100.0

    self.B1_k = a1_k + b1_k*self.acentric_factor + c1_k*self.adimensional_dipole_moment**4.0 + d1_k*self.association_factor
    self.B2_k = a2_k + b2_k*self.acentric_factor + c2_k*self.adimensional_dipole_moment**4.0 + d2_k*self.association_factor
    self.B3_k = a3_k + b3_k*self.acentric_factor + c3_k*self.adimensional_dipole_moment**4.0 + d3_k*self.association_factor
    self.B4_k = a4_k + b4_k*self.acentric_factor + c4_k*self.adimensional_dipole_moment**4.0 + d4_k*self.association_factor
    self.B5_k = a5_k + b5_k*self.acentric_factor + c5_k*self.adimensional_dipole_moment**4.0 + d5_k*self.association_factor
    self.B6_k = a6_k + b6_k*self.acentric_factor + c6_k*self.adimensional_dipole_moment**4.0 + d6_k*self.association_factor
    self.B7_k = a7_k + b7_k*self.acentric_factor + c7_k*self.adimensional_dipole_moment**4.0 + d7_k*self.association_factor

  ### Methods 

  def calculateDynamicViscosity(self, P, T, rho):

    # T. H. Chung, L. L. Lee, K. E. Starling.
    # Applications of kinetic gas theories and multiparameter correlation for prediction of dilute gas viscosity and thermal conductivity.
    # Industrial & Engineering Chemistry Fundamentals, 23, 8-13, 1984.

    # T. H. Chung, M. Ajlan, L. L. Lee, K. E. Starling.
    # Generalized multiparameter correlation for nonpolar and polar fluid transport properties.
    # Industrial & Engineering Chemistry Fundamentals, 27, 671-679, 1988.

    # Auxiliar coefficients

   # rho = self.calculateDensityInternalEnergyFromPressureTemperature( P, T)     ### Call rho calculate in one method of BaseThermodynamics model
    #T   = self.                                                                 ### Call T calculate in one method of BaseThermodynamics model

    v     = self.molecular_weight/rho
    Y     = self.critical_molar_volume/( 6.0*v )
    T_ast = 1.2593*( T/self.critical_temperature )
    Omega = 1.16145*T_ast**(-0.14874)  + 0.52487*np.exp( -0.77320*T_ast ) + 2.16178*np.exp( -2.43787*T_ast )
    G1    = ( 1.0 - 0.5*Y )/((1.0 - Y)**3.0)
    Fc    = 1.0 - 0.2756*self.acentric_factor + 0.059035*self.adimensional_dipole_moment**4.0 + self.association_factor
    
    # Additional auxiliar coefficients

    G2_mu      = ( self.E1_mu*( 1.0 - np.exp( -self.E4_mu*Y ) )/Y + self.E2_mu*G1*np.exp( self.E5_mu*Y ) + self.E3_mu*G1 )/( self.E1_mu*self.E4_mu + self.E2_mu + self.E3_mu )
    mu_ast_ast = (self.E7_mu*Y**2.0)*G2_mu*np.exp( (self.E8_mu + self.E9_mu/T_ast) + self.E10_mu*T_ast**(-2.0))
    mu_ast     = ( np.sqrt( T_ast )*Fc/Omega )*( 1.0/G2_mu + self.E6_mu*Y ) + mu_ast_ast

    # Calculate viscosity
    mu = (1.0e-7)*mu_ast*( ( 36.344*np.sqrt( ( 1.0e3*self.molecular_weight )*self.critical_temperature ) )/ (1.0e6*self.critical_molar_volume)**(2.0/3.0)) 

    return( mu )
  

  def calculateThermalConductivity(self, P, T, rho): 

    # T. H. Chung, L. L. Lee, K. E. Starling.
    # Applications of kinetic gas theories and multiparameter correlation for prediction of dilute gas viscosity and thermal conductivity.
    # Industrial & Engineering Chemistry Fundamentals, 23, 8-13, 1984.

    # T. H. Chung, M. Ajlan, L. L. Lee, K. E. Starling.
    # Generalized multiparameter correlation for nonpolar and polar fluid transport properties.
    # Industrial & Engineering Chemistry Fundamentals, 27, 671-679, 1988.

    # Auxiliar coefficients

    std_bar_c_p =  self.calculateMolarStdCpFromNASApolynomials(T)

    v           = self.molecular_weight/rho
    Y           = self.critical_molar_volume/( 6.0*v )
    T_ast       = 1.2593*( T/self.critical_temperature )
    Omega       = 1.16145*pow( T_ast, -0.14874 ) + 0.52487*np.exp( -0.77320*T_ast ) + 2.16178*np.exp( -2.43787*T_ast )
    G1          = ( 1.0 - 0.5*Y )/pow( 1.0 - Y, 3.0 )
    Fc          = 1.0 - 0.2756*self.acentric_factor + 0.059035*pow( self.adimensional_dipole_moment, 4.0 ) + self.association_factor
    std_bar_c_p = self.calculateMolarStdCpFromNASApolynomials( T )

    # Additional auxiliar coefficients
    mu_0_k  = 40.785e-7*Fc*np.sqrt( 1.0e3*self.molecular_weight*T )/( pow( 1.0e6*self.critical_molar_volume, 2.0/3.0 )*Omega )
    alpha_k = ( std_bar_c_p/self.R_universal - 1.0 ) - 1.5                                                                                   
    beta_k  = 0.7862 - 0.7109*self.acentric_factor + 1.3168*pow( self.acentric_factor, 2.0 )
    gamma_k = 2.0 + 10.5*pow( T/self.critical_temperature, 2.0 )
    Psi_k   = 1.0 + alpha_k*( ( 0.215 + 0.28288*alpha_k - 1.061*beta_k + 0.26665*gamma_k )/( 0.6366 + beta_k*gamma_k + 1.061*alpha_k*beta_k) )
    q_k     = 0.003586*( np.sqrt( self.critical_temperature/self.molecular_weight )/pow( ( 1.0e6*self.critical_molar_volume ), 2.0/3.0 ) )
    G3_k    = ( ( ( self.B1_k/Y )*( 1.0 - np.exp( (-1.0)*self.B4_k*Y ) ) ) + ( self.B2_k*G1*np.exp( self.B5_k*Y ) ) + ( self.B3_k*G1 ) )/( self.B1_k*self.B4_k + self.B2_k + self.B3_k )

    # Calculate thermal conductivity
    kappa = ( 31.2*mu_0_k*Psi_k/self.molecular_weight )*( 1.0/G3_k + self.B6_k*Y ) + q_k*self.B7_k*pow( Y, 2.0 )*np.sqrt( T/self.critical_temperature)*G3_k

    return( kappa )
  

  def calculateMolarStdCpFromNASApolynomials(self, T): 

    std_bar_c_p = 0.0
    if (T>=200.0) and (T<1000.0) :
        std_bar_c_p = self.R_universal*(self.NASA_coefficients[7] + self.NASA_coefficients[8]*T + self.NASA_coefficients[9]*T**2.0 + self.NASA_coefficients[10]*T**3.0 + self.NASA_coefficients[11]* T**4.0)
    elif (T>=1000.0) and (T<6000.0):
        std_bar_c_p = self.R_universal*(self.NASA_coefficients[0] + self.NASA_coefficients[1]*T + self.NASA_coefficients[2]*T**2.0 + self.NASA_coefficients[3]*T**3.0 + self.NASA_coefficients[4]*T**4.0)
    elif (T < 200):
        # Assume constant temperature below T = 200 K	    
        T_min = 200.0	    
        std_bar_c_p = self.R_universal*(self.NASA_coefficients[7] + self.NASA_coefficients[8]*T_min + self.NASA_coefficients[9]*T_min**2.0 + self.NASA_coefficients[10]*T_min**3.0 + self.NASA_coefficients[11]*T_min**4.0)
    else:
        print(f"\nNASA 7-coefficient polynomials for std bar c_p. T = {T} is above 6000 K.\n\n")
        exit()
    return std_bar_c_p
  