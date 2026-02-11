#ifndef POCKETFFT_H
#define POCKETFFT_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

	typedef struct pocketfft_plan_i* pocketfft_plan;
	pocketfft_plan pocketfft_create_plan(int n);
	void pocketfft_destroy_plan(pocketfft_plan plan);
	void pocketfft_execute(pocketfft_plan plan, float* data, float* result);

#ifdef __cplusplus
}
#endif

#endif