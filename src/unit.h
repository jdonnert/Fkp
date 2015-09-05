/* Code Units */
void set_units();

struct units {			/* For unit Conversion. */
	double Length;
	double Mass;
	double Vel;
    double Time;
    double Energy;
} Unit, Comv2phys;

double temperature_cgs(float);
double density_cgs(float);
double number_density_cgs(float);
double proton_number_density_cgs(float);
double thermal_energy_density_cgs(float, float);
