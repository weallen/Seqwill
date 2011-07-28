#include "analysis/Random.h"

gsl_rng* init_rng() {
	const gsl_rng_type* T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	return gsl_rng_alloc(T);
}
