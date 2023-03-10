//============================================================================
// Name        : Density_HV_and_SG_Cal.cpp
// Author      : Abdurrahman Nurhakim
// Version     : 1.0
// Copyright   : Your copyright notice
// Description : Calculate Density and Spesific Gravity in C++
//============================================================================

#include <stdio.h>
#include "Gross.h" //this library comes from Aga8 standard library, you can find this library at https://github.com/usnistgov/AGA8/blob/master/AGA8CODE/C/Gross.h
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <Detail.h>  //this library comes from Aga8 standard library, you can find this library at https://github.com/usnistgov/AGA8/blob/master/AGA8CODE/C/Detail.h
#include <Density_HV_and_SG_Cal.h>

Density_HV_and_SG_Cal::Density_HV_and_SG_Cal() {

}

float Density_HV_and_SG_Cal::to_kg_per_m3(float mol_per_m3, float kg_per_mol) {
	float unit_kg_per_m3;
	unit_kg_per_m3 = mol_per_m3 * kg_per_mol;
	return unit_kg_per_m3;
}

void Density_HV_and_SG_Cal::molar_decimal(struct Molar_Gas Gas_Mol_conf,
		double *_molar) {
	_molar[0] = (double) Gas_Mol_conf._Methane / 100.0;
	_molar[1] = (double) Gas_Mol_conf._Nitrogen / 100.0;
	_molar[2] = (double) Gas_Mol_conf._Carbon_Dioxide / 100.0;
	_molar[3] = (double) Gas_Mol_conf._Ethane / 100.0;
	_molar[4] = (double) Gas_Mol_conf._Propane / 100.0;
	_molar[5] = (double) Gas_Mol_conf._I_Butane / 100.0;

	_molar[6] = (double) Gas_Mol_conf._N_Butane / 100.0;
	_molar[7] = (double) Gas_Mol_conf._I_Pentane / 100.0;
	_molar[8] = (double) Gas_Mol_conf._N_Pentane / 100.0;
	_molar[9] = (double) Gas_Mol_conf._N_Hexane / 100.0;

	_molar[10] = (double) Gas_Mol_conf._N_Heptane / 100.0;
	_molar[11] = (double) Gas_Mol_conf._N_Octane / 100.0;
	_molar[12] = (double) Gas_Mol_conf._N_Nonane / 100.0;
	_molar[13] = (double) Gas_Mol_conf._N_Decane / 100.0;

	_molar[14] = (double) Gas_Mol_conf._Hydrogen / 100.0;
	_molar[15] = (double) Gas_Mol_conf._Oxygen / 100.0;
	_molar[16] = (double) Gas_Mol_conf._Carbon_Monoxide / 100.0;
	_molar[17] = (double) Gas_Mol_conf._Water / 100.0;
	_molar[18] = (double) Gas_Mol_conf._Hydrogen_Sulfide / 100.0;

	_molar[19] = (double) Gas_Mol_conf._Hellium / 100.0;
	_molar[20] = (double) Gas_Mol_conf._Argon / 100.0;
}

void Density_HV_and_SG_Cal::DensityCalculate(struct Molar_Gas Gas_Mol_conf,
		struct PTD_Parameter data_in, float *output) {
	struct out_detail_aga8 aga8_detail_result_;
	double total, density;
	float _P1 = data_in.Pressure * 100000.0;	//Upstream Pressure (Pa G)
	float _T = data_in.Temperature + 273.15;	//Downstream Temperature (K)	
	double *_molar = (double*) malloc(sizeof(double) * 21); //dinamic memory allocation
	double _K = (double) _T;
	double _kPa = (double) (_P1) / 1000.0;
	this->molar_decimal(Gas_Mol_conf, _molar);
	aga8_detail_result_ = calculate_detail_aga8(_molar, _K, _kPa);
	density = this->to_kg_per_m3(aga8_detail_result_.Molar_density,
			aga8_detail_result_.Molar_mass); ///density
	free(_molar); //free
	*output = density;
}

struct _output_Gross_Aga8 Density_HV_and_SG_Cal::calculate_Gross_Aga8(
		struct Molar_Gas Gas_Mol_conf, double temp, double pressure,
		double base_temp, double base_pressure) {
	double temp_in_K_h = 298.15;
	double OFFSET_BASE_TEMPERATURE;
	double temp_in_K_d = base_temp + 273.15;
	double temp_in_K = temp + 273.15;
	double press_in_kPa_d = base_pressure * 100;
	double press_in_kPa = (pressure) * 100;
	_output_Gross_Aga8 buff_output_Gross_Aga8;
	double *_x = (double*) malloc(sizeof(double) * 21);
	this->molar_decimal(Gas_Mol_conf, _x);
	SetupGross();
	std::vector<double> x(_x, _x + 21), xGrs(4, 0);
	x.insert(x.begin(), 0.0);
	double mm = 0;
	MolarMassGross(x, mm);
	int ierr = 0;
	std::string herr;
	double T = temp_in_K, P = press_in_kPa, Td = temp_in_K_d, Th = temp_in_K_h,
			Pd = press_in_kPa_d, D = 6.35826, pp = -1, Hv = -1, Hv2 = -1, Gr,
			HN, HCH, Z = -1;

	printf("\n\nmine-Inputs-----\n");
	printf("Temperature [K]:                    300.0000000000000 != %0.16g\n",
			T);
	printf("Pressure [kPa]:                     10000.00000000000 != %0.16g\n",
			P);
	printf("Th [K]:                             298.1500000000000 != %0.16g\n",
			Th);
	printf("Td [K]:                             273.1500000000000 != %0.16g\n",
			Td);
	printf("Pd [kPa]:                           101.3250000000000 != %0.16g\n",
			Pd);

	GrossHv(x, xGrs, HN, HCH);
	DensityGross(T, P, xGrs, HCH, D, ierr, herr);
	PressureGross(T, D, xGrs, HCH, pp, Z, ierr, herr);

	printf("Outputs-----\n");
	printf("Molar mass [g/mol]:                 20.54333051000000 != %0.16g\n",
			mm);
	printf("Molar density [mol/l]:              5.117641317088482 != %0.16g\n",
			D);
	printf("Pressure [kPa]:                     9999.999999999998 != %0.16g\n",
			P);
	printf("Compressibility factor:             0.7833795701012788 != %0.16g\n",
			Z);

	GrossInputs(Td, Pd, x, xGrs, Gr, HN, HCH, ierr, herr);
	GrossHv(x, xGrs, HN, HCH);
	DensityGross(Td, Pd, xGrs, HCH, D, ierr, herr);
	Hv = HN * D;

	printf("Molar density [mol/l]				: 5.117641317088482 != %0.16g\n", D);
	printf(
			"Relative density (Spesific Gravity)	: 0.7112387718599272 != %0.16g\n",
			Gr);
	printf("Molar heating value (kJ/mol, 25 C)	: 924.3591780000000 != %0.16g\n",
			HN);
	printf("HCH (kJ/mol, 25 C)					: 1004.738236956522 != %0.16g\n", HCH);
	printf("Equivalent hydrocarbon fraction		: 0.9199999999999999 != %0.16g\n",
			xGrs[1]);
	printf("nitrogen mole fraction				: 2.000000000000000E-02 != %0.16g\n",
			xGrs[2]);
	printf("CO2 mole fraction					: 6.000000000000000E-02 != %0.16g\n",
			xGrs[3]);
	printf("Volumetric heating value at Td,Pd	:  41.37676728724603 != %0.16g\n",
			Hv * 26.86274632);

	xGrs[2] = x[2];
	xGrs[3] = x[3];
	GrossMethod2(Th, Td, Pd, xGrs, Gr, Hv2, mm, HCH, HN, ierr, herr);
	DensityGross(T, P, xGrs, HCH, D, ierr, herr);

	printf(
			"Volumetric heating value at Td,Pd 	:  41.37676728724603 != %0.16g\n",
			Hv2 * 26.86274632);
	printf("Gross method 2-----\n");
	printf("Molar density [mol/l]				: 5.197833636353455 != %0.16g\n", D);

	xGrs[2] = x[2];
	xGrs[3] = x[3];
	GrossMethod1(Th, Td, Pd, xGrs, Gr, Hv, mm, HCH, HN, ierr, herr);
	DensityGross(T, P, xGrs, HCH, D, ierr, herr);

	printf("Gross method 1-----\n");
	printf(
			"Molar density [mol/l]:              5.144374668159809 != %0.16g\n\n\n",
			D);

	buff_output_Gross_Aga8.Molar_mass = mm;
	buff_output_Gross_Aga8.Molar_density = D;
	buff_output_Gross_Aga8.Gr = Gr; ////Spesific Gravity
	buff_output_Gross_Aga8.HN = HN;
	buff_output_Gross_Aga8.HCH = HCH;
	buff_output_Gross_Aga8.xGrs[1] = xGrs[1];
	buff_output_Gross_Aga8.xGrs[2] = xGrs[2];
	buff_output_Gross_Aga8.xGrs[3] = xGrs[3];
	buff_output_Gross_Aga8.Hv = Hv; ////Heating Value Method 1
	buff_output_Gross_Aga8.Hv2 = Hv2; ////Heating Value Method 2

	///Constant Value Based on -> https://www.unitsconverters.com/en/Mj/M3-To-Btu/Ft3/Utu-4795-4789
	buff_output_Gross_Aga8.BTU_SCF = Hv2 * 26.86274632; ////Heating Value Method 2 in BTU/SCF Unit
	free(_x);
	return buff_output_Gross_Aga8;
}



/* FAIL !!!
 void Density_HV_and_SG_Cal::spesific_gravity(struct Molar_Gas Gas_Mol_conf, struct ABBBS_Parameter_SG str_SG_, int type_input, float *output){
 struct PTD_Parameter str_;
 switch(type_input){
 case _MANUAL_INPUT:
 str_.Pressure = str_SG_.Base_Pressure;
 str_.Temperature = str_SG_.Base_Temperature;
 this->DensityCalculate(Gas_Mol_conf, str_, &str_SG_.Base_Density);
 str_SG_.Spesific_Gravity = str_SG_.Base_Density / str_SG_.Air_Density;  //https://www.sciencedirect.com/topics/materials-science/density-specific-gravity
 *output = str_SG_.Spesific_Gravity;
 break;
 case _LIBRARY_INPUT:
 str_SG_.Air_Density = constant_air_density;
 str_.Pressure = constant_base_pressure;
 str_.Temperature = constant_base_temperature;
 this->DensityCalculate(Gas_Mol_conf, str_, &str_SG_.Base_Density);
 str_SG_.Spesific_Gravity = str_SG_.Base_Density / str_SG_.Air_Density; //https://www.sciencedirect.com/topics/materials-science/density-specific-gravity
 *output = str_SG_.Spesific_Gravity;
 break;
 }
 }
 */

/***
 ////HOW TO USE IT

 ////It also can be used in Arduino IDE

 ////example code (eclipse) :
 int main (){
 struct _output_Gross_Aga8 Output_Parameter;
 struct PTD_Parameter _data_in;
 struct Molar_Gas Gas_Mol_conf;
 Density_HV_and_SG_Cal ClassLibrary;
 double _pressure, _temperature; //Input Pressure and Temperature
 double _base_pressure, _base_temperature; ///Base Parameter
 float molar_[21] = {77.824, 2, 6, 8, 3, 0.15, 0.3,
 0.05, 0.165, 0.215, 0.088, 0.024, 0.015,
 0.009, 0.4, 0.5, 0.2, 0.01, 0.0025, 0.7, 0.1 };
 float Density;
 
 while(1) { ///LOOP FOREVER (you also can use --> for(;;))
 pressure = 22.00; //in bar
 temperature = 22.00; //in celcius
 base_temperature = 15; //in celcius
 base_pressure = 1.01325; //in bar (athmosphere pressure)

 Gas_Mol_conf._Methane = molar_[0];
 Gas_Mol_conf._Nitrogen = molar_[1];
 Gas_Mol_conf._Carbon_Dioxide = molar_[2];
 Gas_Mol_conf._Ethane = molar_[3];
 Gas_Mol_conf._Propane = molar_[4];
 Gas_Mol_conf._I_Butane = molar_[5];
 Gas_Mol_conf._N_Butane = molar_[6];
 Gas_Mol_conf._I_Pentane = molar_[7];
 Gas_Mol_conf._N_Pentane = molar_[8];
 Gas_Mol_conf._N_Hexane = molar_[9];
 Gas_Mol_conf._N_Heptane = molar_[10];
 Gas_Mol_conf._N_Octane = molar_[11];
 Gas_Mol_conf._N_Nonane = molar_[12];
 Gas_Mol_conf._N_Decane = molar_[13];
 Gas_Mol_conf._Hydrogen = molar_[14];
 Gas_Mol_conf._Oxygen = molar_[15];
 Gas_Mol_conf._Carbon_Monoxide = molar_[16];
 Gas_Mol_conf._Water = molar_[17];
 Gas_Mol_conf._Hydrogen_Sulfide = molar_[18];
 Gas_Mol_conf._Hellium = molar_[19];
 Gas_Mol_conf._Argon = molar_[20];

 _pressure = 22.00; //in bar
 _temperature = 22.00; //in celcius
 _base_temperature = 15; //in celcius
 _base_pressure = 1.01325; //in bar (athmosphere pressure)

 ClassLibrary.DensityCalculate(Gas_Mol_conf, _data_in, &Density);
 Output_Parameter = ClassLibrary.calculate_Gross_Aga8(Gas_Mol_conf,
 _temperature, _pressure, _base_temperature, _base_pressure); //Gross.h
 printf("Heating Value: %f\n\n", Output_Parameter.BTU_SCF);
 printf("Spesific Grafity: %f\n\n", Output_Parameter.Gr);
 printf("Density: %f\n\n", Density);

 }
 
 }

 ***/

