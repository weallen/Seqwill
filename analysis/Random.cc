#include "analysis/Random.h"

void Rng::Init() {
	const gsl_rng_type* T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng_ = gsl_rng_alloc(T);
}

