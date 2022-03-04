**************************************************************************
* This file presents a model script for Breakdown simulations
*	Author: Sinuo Zhang
* 	Last Modified:Fri 04 Mar 2022 04:07:38 PM CET
*	Description:
*	The breakdown simulation of the structures with floating
*		implants or using materials with large band-
*		gaps, requires well-chosen methodes for the breakdown
*		simulaiton.
*	Surface damge model: 	Perugia (2021)
*	Bulk damage model:		Hamburg Penta Trap Model (2018)
*	Parameters need to be specified in the workbench:
*		> BS		:	Bias Voltage (V) -> positive float
*		> TIDMrad	:	The TID dose (Mrad)
*		> Neq		:	Neutron equivalent fluence (cm^-2)
*		> Model		: 	Avalanche generation model -> name
*		> Temp		:	Temperature of the device (K)
*		> transport :	Transportation model -> see code
**************************************************************************

* define the maximum ramping current for the break criteria, used in "Math" section
#define _Imax_ 2.0e-9
* The resistor which attatched to the electrode as the bias resistor. This is essential for the "Resistor method" of the breakdown simulation.
* It is also used in the "transient method" breakdown simulation
#define _Resistor 1e+4
* In this case we use the transient simulation to ramp up the bias voltage, since this method is more robust for structures with floating implants
* 	Thus we set an end time for ramping, which can be e.g. the bias voltage
*	It should be noticed that the ramping time(speed) can affect the result, and tests should be done before increasing the simulation speed
#define _tmax_ @BS@

* Define the electrodes
	{name= "NW"		Voltage=0.0 }
	#if @NRCont@ == 1
		{name= "NR"	Voltage=0.0 }
	#endif
	{name= "BS"	 	Voltage=0.0 Voltage= (0.0 at 0.0, -@BS@ at _tmax_) Resistor= _Resistor}
}

File {
	Grid		=	"@tdr@"
	Plot		=	"@tdrdat@"
	Parameter	=	"@parameter@"
	Current		=	"@plot@"
}

Physics {
	*Using Fermi Statistics
	Fermi
	* The default transportation model is the "Drift Diffusion" model, which doesn't need to be specified
	* The "Hydrodynamic" model takes the charge carrier temperature into account. The results can be more accurate for breakdown simulations, since high current can occur and the temperature is an important influencing factor.
	#if "@transport@" == "hd"
		*The specified carrier temperature will be included in the model
		Hydrodynamic (eTemperature hTemperature)
	#endif
	*The band gap model: Slotboom (Temperature dependent band gap size)
	*e.g. the Slotboom is also available for affinity model
	*The bandgap model is the sum of E0 and Ef parts, where the E0 is determined
	* by such models, and the Ef is (by default switched on for Fermi statisics)
	* an option to correct the band gap obtained via Maxwell-Boltzmann statistics
	* (an oftenly used assumption for parameter extraction from experimental data, 
	* larger error for high doping regim).
	*This correction is more desirable for the simulation without using
	* Fermi statistics
	*The Ef can be switched off by specifying "NoFermi"
	*Sometimes for III-IV materials the "fermi correction" is too large
	Temperature= @Temp@
	EffectiveIntrinsicDensity(BandGapNarrowing(OldSlotboom))
	*Mobility models
	Mobility (
		*Mobility degradation due to impurity scattering --> doping dependent
		*There exist sveral models and options, diff. material has diff. default
		* model, for silicon is "Masetti", but (the default model) can also be 
		* set in the parameter file.
		*If more models are used, they are combined by the matthiessen's rule
		* -> reciprocal sum.
		DopingDependence(Masetti)
		*e-Velocity is no more proportional to the electric field at high field
		* it saturates at a certain value.
		*In orther words, the mobility is not constant any more.
		*The default is Canali model
		HighFieldSaturation
		*The "ConwellWeisskopf" is a model based on Conwell-Weisskopf screening
		CarrierCarrierScattering(ConwellWeisskopf)
	)
	*Generation and recombination of charge carriers
	Recombination (
		*Shockley-Read-Hall recombination model (G/R through traps in bandgap)
		SRH(
			*SRH-lifetime depends on doping concentration, which is modeled
			* with the Scharfetter relation
			*The lifetime can also be specified via 
			DopingDependence
			*Temperature dependent lifetime
			*There are power and exponential law of the dependence
			*Temperature behaviour depends strongly on the nature of 
			* recombination centres --> no universial law
			TempDependence
			*Electric field driven trap-assisted tunneling with Hurkx model
			*Describes local field-dependence of the SRH lifetime
			*For converging problem: specify "NoSRHperPotential" in math section
			*The tunneling mass of the model can be specified in parameter file
			Tunneling(Hurkx)
			ElectricField(Lifetime=Hurkx DensityCorrection= None)
		)
		*The auger recombination is important at high carrier densities
		*The Auger coefficient is by default positive --> pure recombination
		*The Generation can be specified by "WithGeneration"
		# Auger(WithGeneration)
		*Avalanche effect (impact ionisation)
		*Default is the "vanOverstreaten" model
		*Can also be specified as driving force model
		*Default driving force is the "GradQuasiFermi" gradient of the
		* quasi-Fermi potential

		Avalanche(@Model@ Eparallel)

		*Band to band recombination according to the Hurkx model
		Band2Band(Hurkx)
	)	
}

* ----- Begin of the TID surface damage section ------
# The template can be found in the appendix of the article
*Dose of the TID in unit of Mrad
#set TID 		[format %.2f @TIDMrad@] 
* Initial total conctration of the oxide charge before irradiation
* since we don't know it for LF, we use a dummy value
#set Qoxpre 	[format %.2e 1.0e+10]		
* Initial total acceptor concentration at the interface before irradiation
* since we don't know it for LF, we use a dummy value
#set Nitaccpre 	[format %.2e 1.0e+9]	
* Initial total donor concentration at the interface before irradiation
* since we don't know it for LF, we use a dummy value
#set Nitdonpre 	[format %.2e 1.0e+9]		

* Case selection: before TID
#if "@TIDMrad@" == 0						
	# due to TID introduced total damage concentration (for traps indicated with "N" [cm^-2]) are set to 0
	#set DeltaQox 		[format %.2e 0]
	#set DeltaNitacc 	[format %.2e 0]
	#set DeltaNitdon 	[format %.2e 0]
	# the energy concentration of the traps labelled with "D" [eV^-1 cm^-2].
	# As we see also later by specifying the traps, this concentration is used for the uniform distribution of the trap energies
	#	therefore this value represents the height of the distribution.
	# Converting the "D" to "N" requires an integration of "D" over the energy. since it's a uniform distribution, a simple
	# 	conversion function is: N=D*width(eV) and D=N/width(eV).
	# Hence the definition and conversion beneath.
	* [SNOW: for Dit_acc, the original "/0.3" is replaced by "0.56", s.u.]
	#set Dit_acc 		[format %.2e @<Nitaccpre/0.56>@]		
	#set Dit_don		[format %.2e @<Nitdonpre/0.3>@] 
	#set Qox 			[format %.2e @Qoxpre@]
#else						
* Case seletion: after TID
	# due to TID introduced trap ("N") and oxide charge concentration are added according to fitting the data from measurements
	#set DeltaQox 		[format %.2e [expr @<3.74E+11 + 6.20E+10 * log(TID)>@]] 
	#set DeltaNitacc 	[format %.2e [expr @<6.35E+11 + 1.50e+11 * log(TID)>@]] 
	#set DeltaNitdon 	[format %.2e [expr @<1.07E+12 + 2.90e+11 * log(TID)>@]]

	# The total damage parameters (the "EsA" is probably a typo, SNOW: I suggest this should be 0.56)
	#set Dit_acc [format %.2e @<(DeltaNitacc+Nitaccpre)/0.56>@] 
	#set Dit_don [format %.2e @<(DeltaNitdon+Nitdonpre)/0.3>@] 
	#set Qox [format %.2e @<Qoxpre+DeltaQox>@]
#endif

	Physics (MaterialInterface= "Oxide/Silicon") {
	Traps(
		(FixedCharge Conc=@Qox@)
		(Acceptor Conc=@Dit_acc@ Uniform EnergyMid=0.84 EnergySig=0.56 fromValBand eXsection=1e-16 hXsection=1e-15 Add2TotalDoping)
		(Donor Conc=@Dit_don@ Uniform EnergyMid=0.60 EnergySig=0.30 fromValBand eXsection=1e-15 hXsection=1e-16 Add2TotalDoping)
	)
	}
*------ End of the Surface damage section ------

*------ Begin of the Bulk damage section ------
	Physics (Material="Silicon") {
	Traps (
		* ===E30K===
		( Donor Level fromCondBand Conc=@<Neq*0.0497>@ EnergyMid=0.1 eXsection=2.3E-14 hXsection=2.92E-16 )
		* ====V3======         		
		( Acceptor Level fromCondBand Conc=@<Neq*0.6447>@ EnergyMid=0.458 eXsection=2.551E-14 hXsection=1.511E-13 )
		* =====Ip======
		( Acceptor Level fromCondBand Conc=@<Neq*0.4335>@ EnergyMid=0.545 eXsection=4.478E-15 hXsection=6.709E-15 )
		* ======H220=====
		( Donor Level fromValBand Conc=@<Neq*0.5978>@ EnergyMid=0.48 eXsection=4.166E-15 hXsection=1.965E-16)
		* ======CiOi======
		( Donor Level fromValBand Conc=@<Neq*0.3780>@ EnergyMid=0.36 eXsection=3.23E-17 hXsection=2.036E-14)
	)
	}
*------ End of the Bulk damage section -----

Plot {
	*- doping profiles
	Doping 
	# DonorConcentration AcceptorConcentration
	*- charge field potential and potential energy
	SpaceCharge
	ElectricField/Vector Potential
	# BandGap EffectiveBandGap BandgapNarrowing ElectronAffinity
	ConductionBandEnergy ValenceBandEnergy
	*- carrier densities
	# EffectiveIntrinsicDensity IntrinsicDensity
	eDensity hDensity
	eQuasiFermiEnergy hQuasiFermiEnergy
	# HeavyIonChargeDensity
	*- Temperatures
	eTemperature hTemperature Temperature
	*- current and current components
	eGradQuasiFermi/Vector hGradQuasiFermi/Vector
	eMobility hMobility eVelocity hVelocity
	Current/Vector eCurrent/Vector hCurrent/Vector
	eDriftVelocity/Vector hDriftvelocity/Vector
	*- SRH & interface traps
	*	"tSRHRecombination" recombination rate at the defect level
	SRHRecombination
	tSRHRecombination
	eLifetime hLifetime
	*- Band2Band tunneling
	eBand2BandGeneration hBand2BandGeneration Band2BandGeneration
	*- Avalanche effect
	eAvalancheGeneration hAvalancheGeneration AvalancheGeneration
	eIonIntegral HIonIntegral MeanIonIntegral
	eAlphaAvalanche hAlphaAvalanche
}

Math {
	* Use Extrapolation in numerical simulation to provide a approximate prediction of the next value
	Extrapolate
	* Control on relative and absolute errors, This enables the relative error control
	RelErrControl
	* Force the convergence criteria to be: both RHS and Error are smaller than the specified minimum reguired value
	*	In general, a more strict convergence criteria helps to improve the overall convergence performance
	*	This is a strict criteria, nevertheless the exact performance need to be tested
	# RhsAndUpdateConvergence
	* force additional checking on RHS norm
	CheckRhsAfterUpdate
	* Set the minimum value of RHS, below which the calculation is considered as converged
	RhsMin= 1e-10
	* relative error = 10^(-Digits)
	Digits= 8
	* use full derivatives in Newton method
	Derivatives
	* use avalanche derivatives
	AvalDerivatives
	* relavtive error
	ErrRef(electron)=1.e8
	ErrRef(hole)=1.e8
	* Precision
	ExtendedPrecision
	* Set the breakdown criterion, i.e break when the current above a limit value
	#if @NRCont@ == 0
	BreakCriteria{ Current(Contact= "NW" AbsVal= _Imax_) } 
	#elif @NRCont@ == 1
	BreakCriteria{ Current(Contact= "NR" AbsVal= _Imax_) }
	#endif
	* Transient method, In this case we use the "BE" (Backwards Euler method), which is the more rubost method in sentaurus TCAD. It's accuracy may be worse than the default method. Nevertheless, a less accurate result is better than nothing comming out :)
	Transient= BE
	*numerical parameter for space-charge regions
	#eDrForceRefDens= 1e10
	#hDrForceRefDens= 1e10
	* Maximum number of iterations at each step
	Iterations= 15	* Too many iterations may benefit the convergence of each step. However it may not be time efficient.
					* 	Therefore it's recommended to use less iterations per step, and let the solver use smaller step size
	NotDamped= 100  * This is larger than iterations, meaning better not to trigger damping.
					*	It is recommended to do so by experts and official examples.
					*	Damping is perhaps beneficial in special cases
	* Solver of the linear system
	* The Default solver "Super" is the most robust solver, nevertheless it doesn't support multi-threading -> slow but good
	*	The "ParDiso" supports multi-threading, but it's not so robust. However, when it converges it's still okay
	*	The robutness of "ILS" lies between "Super" and "ParDiso" (according to experts), it also support multi-threading
	Method= ParDiso
	*display simulation time in "human" units
	Wallclock
	*display max.error information
	#CNormPrint
	* to avoid convergence problem when simulating defect-assisted tunneling
	NoSRHperPotential
	RecBoxIntegr
	* Set multi-threading (not working for default linear solver)
	NumberOfThreads=10
}

Solve {
		* Provide an initial value of the simulation. Since this has only one step, a large iteration can be used to ensure the convergence
		Coupled(Iterations=100){Poisson}
		* Provide the initial value for coupled Poission (potential and field) & the charge carriers. For some simulation models, the temperature should also be specified
		#if "@transport@" == "dd"
			Coupled{Poisson Electron Hole}
		#elif "@transport@" == "hd"
			Coupled{Poisson Electron Hole eTemperature hTemperature}
		#endif
		
	#if @BS@ != 0
		* Define the property of the Transient simulation
		* 	The simulaiton speed, convergence properties can be influenced by the values.
		Transient (
			InitialTime= 0.0 FinalTime= _tmax_
			InitialStep= 1e-3
			MaxStep= 2.0	MinStep= 1e-8
			Increment=1.3	Decrement=2.0
		)
		{
		* Solve the equations
		#if "@transport@" == "dd"
		Coupled{ Poisson Electron Hole }
		#elif "@transport@" == "hd"
		Coupled{Poisson Electron Hole eTemperature hTemperature}
		#endif
		* Plot the intermediate tdr of the structure for inspections
		Plot ( FilePrefix = "n@node@_BV" Time = (Range = (0 _tmax_) Intervals=9) NoOverwrite)
		}
		#endif
}	
