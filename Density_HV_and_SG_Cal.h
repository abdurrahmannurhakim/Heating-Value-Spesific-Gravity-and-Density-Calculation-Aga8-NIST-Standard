//============================================================================
// Name        : Density_HV_and_SG_Cal.h
// Author      : Abdurrahman Nurhakim
// Version     : 1.0
// Copyright   : Your copyright notice
// Description : Calculate Density and Spesific Gravity in C++
//============================================================================

#include <stdlib.h>
#include "Gross.h" ////this library comes from Aga8 standard library, you can find this library at https://github.com/usnistgov/AGA8/blob/master/AGA8CODE/C/Gross.h
#include <stdint.h>
#include <Detail.h>  //this library comes from Aga8 standard library, you can find this library at https://github.com/usnistgov/AGA8/blob/master/AGA8CODE/C/Detail.h
#define constant_air_density 1.293 //air density from NASA: https://www.earthdata.nasa.gov/topics/atmosphere/atmospheric-pressure/air-mass-density
#define constant_base_pressure 1.01325 //atmosphere pressure in bar
#define constant_base_temperature 15 //in celcius

struct Molar_Gas{
	float _Methane; 
	float _Nitrogen; 
	float _Carbon_Dioxide; 
	float _Ethane; 
	float _Propane; 
	float _Water; 
	float _Hydrogen_Sulfide; 
	float _Hydrogen; 
	float _Carbon_Monoxide; 
	float _Oxygen; 
	float _I_Butane; 
	float _N_Butane; 
	float _I_Pentane; 
	float _N_Pentane; 
	float _N_Hexane; 
	float _N_Heptane; 
	float _N_Octane; 
	float _N_Nonane; 
	float _N_Decane; 
	float _Hellium; 
	float _Argon; 
};

struct _output_Gross_Aga8 {
	double Molar_mass;
	double Molar_density;
	double Gr; ///Relative density or Spesific Gravity
	double HN; ///Molar heating value (kJ/mol, 25 C)
	double HCH; ///(kJ/mol, 25 C):
	double xGrs[4]; ///quivalent hydrocarbon fraction
	double Hv; ///Volumetric heating value at Td,Pd
	double Hv2; ///Volumetric heating value at Td,Pd Method2
	double BTU_SCF; ///Volumetric heating value at Td,Pd
};

struct PTD_Parameter {
	float Pressure; 
	float Temperature;
	float Density;
};

struct ABBBS_Parameter_SG {
	float Air_Density;
	float Base_Pressure;
	float Base_Temperature;
	float Base_Density;
	float Spesific_Gravity;
};

typedef enum{
 _MANUAL_INPUT = 1,
 _LIBRARY_INPUT
};

class Density_HV_and_SG_Cal{
public:
	Density_HV_and_SG_Cal(); // Constructors
	void DensityCalculate(struct Molar_Gas Gas_Mol_conf, struct PTD_Parameter data_in, float *output);
	//void spesific_gravity(struct Molar_Gas Gas_Mol_conf, struct ABBBS_Parameter_SG data_in, int type_input, float *output);
	void molar_decimal(struct Molar_Gas Gas_Mol_conf, double *out_molar);
	float to_kg_per_m3(float mol_per_m3, float kg_per_mol);
	struct _output_Gross_Aga8 calculate_Gross_Aga8(struct Molar_Gas Gas_Mol_conf, double temp_in_K,
		double press_in_kPa, double base_temp, double base_pressure);  ///Gross.h (SG and Heating Value)
};
