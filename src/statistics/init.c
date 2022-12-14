#include "common.h"
#include "domain.h"
#include "statistics.h"
#include "arrays/statistics.h"


#if NDIMS == 2

/**
 * @brief allocate and initialise statistics_t
 * @param[in] domain : information about domain decomposition and size
 * @return            : pointer to the allocated and initialised statistics_t
 */
statistics_t *statistics_init(const domain_t *domain){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  /* ! structure is allocated ! 1 ! */
  statistics_t *statistics = common_calloc(1, sizeof(statistics_t));
  /* ! arrays are allocated ! 6 ! */
  statistics->ux1   = common_calloc(  UX1_SIZE_0 *   UX1_SIZE_1, sizeof(double));
  statistics->ux2   = common_calloc(  UX2_SIZE_0 *   UX2_SIZE_1, sizeof(double));
  statistics->uy1   = common_calloc(  UY1_SIZE_0 *   UY1_SIZE_1, sizeof(double));
  statistics->uy2   = common_calloc(  UY2_SIZE_0 *   UY2_SIZE_1, sizeof(double));
  statistics->temp1 = common_calloc(TEMP1_SIZE_0 * TEMP1_SIZE_1, sizeof(double));
  statistics->temp2 = common_calloc(TEMP2_SIZE_0 * TEMP2_SIZE_1, sizeof(double));
  statistics->vof1  = common_calloc( VOF1_SIZE_0 *  VOF1_SIZE_1, sizeof(double));
  statistics->vof2  = common_calloc( VOF2_SIZE_0 *  VOF2_SIZE_1, sizeof(double));
  return statistics;
}

#else // NDIMS == 3

/**
 * @brief allocate and initialise statistics_t
 * @param[in] domain : information about domain decomposition and size
 * @return            : pointer to the allocated and initialised statistics_t
 */
statistics_t *statistics_init(const domain_t *domain){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! structure is allocated ! 1 ! */
  statistics_t *statistics = common_calloc(1, sizeof(statistics_t));
  /* ! arrays are allocated ! 8 ! */
  statistics->ux1   = common_calloc(  UX1_SIZE_0 *   UX1_SIZE_1 *   UX1_SIZE_2, sizeof(double));
  statistics->ux2   = common_calloc(  UX2_SIZE_0 *   UX2_SIZE_1 *   UX2_SIZE_2, sizeof(double));
  statistics->uy1   = common_calloc(  UY1_SIZE_0 *   UY1_SIZE_1 *   UY1_SIZE_2, sizeof(double));
  statistics->uy2   = common_calloc(  UY2_SIZE_0 *   UY2_SIZE_1 *   UY2_SIZE_2, sizeof(double));
  statistics->uz1   = common_calloc(  UZ1_SIZE_0 *   UZ1_SIZE_1 *   UZ1_SIZE_2, sizeof(double));
  statistics->uz2   = common_calloc(  UZ2_SIZE_0 *   UZ2_SIZE_1 *   UZ2_SIZE_2, sizeof(double));
  statistics->temp1 = common_calloc(TEMP1_SIZE_0 * TEMP1_SIZE_1 * TEMP1_SIZE_2, sizeof(double));
  statistics->temp2 = common_calloc(TEMP2_SIZE_0 * TEMP2_SIZE_1 * TEMP2_SIZE_2, sizeof(double));
  statistics->vof1  = common_calloc( VOF1_SIZE_0 *  VOF1_SIZE_1 *  VOF1_SIZE_2, sizeof(double));
  statistics->vof2  = common_calloc( VOF2_SIZE_0 *  VOF2_SIZE_1 *  VOF2_SIZE_2, sizeof(double));
  return statistics;
}

#endif // NDIMS
