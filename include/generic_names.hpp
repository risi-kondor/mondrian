/* 
 	Note about function aliasing in C++11: 
 	- 'const auto&' method only works for non-template non-overloaded functions. 
 	- For overloaded functions, use: 
 		const auto& new_fn_name = static_cast<OVERLOADED_FN_TYPE>(old_fn_name);
	- Use "auto g = std::mem_fn(&A::f);" to alias a member function 'f' in class 'A'
	- In C++14, one can alias templated functions as
		template<typename T> constexpr auto alias_to_old = old_function<T>;
*/

/*
 * Headers
 */ 

/*
 * Constants
 */

const auto& BlasColMajor = CblasColMajor;
// const auto& BlasRowMajor = CblasRowMajor;
const auto& BlasNoTrans = CblasNoTrans;
// const auto& BlasTrans = CblasTrans;
// const auto& BlasConjTrans = CblasConjTrans;
// const auto& BlasUpper = CblasUpper;
// const auto& BlasLower = CblasLower;
// const auto& BlasNonUnit = CblasNonUnit;
// const auto& BlasUnit = CblasUnit;
// const auto& BlasLeft = CblasLeft;
// const auto& BlasRight = CblasRight;

/*
 * BLAS Level 1 Routines (vector-vector)
 */

// const auto& blas_drotg = cblas_drotg;
// const auto& blas_drotmg = cblas_drotmg;
// const auto& blas_drot = cblas_drot;
// const auto& blas_drotm = cblas_drotm;
// const auto& blas_dswap = cblas_dswap;
const auto& blas_dscal = cblas_dscal;
// const auto& blas_dcopy = cblas_dcopy;
const auto& blas_daxpy = cblas_daxpy;
const auto& blas_ddot  = cblas_ddot;
// const auto& blas_ddotu = cblas_ddotu;
// const auto& blas_ddotc = cblas_ddotc;
// const auto& blas_dsdot = cblas_dsdot;
const auto& blas_dnrm2 = cblas_dnrm2;
const auto& blas_dasum = cblas_dasum;
const auto& blas_idamax= cblas_idamax;

// extra
const auto& blas_idamin = cblas_idamin;

/*
 * BLAS Level 2 Routines (matrix-vector)
 */

const auto& blas_dgemv = cblas_dgemv;
// const auto& blas_dgbmv = cblas_dgbmv;
// const auto& blas_dhemv = cblas_dhemv;
// const auto& blas_dhbmv = cblas_dhbmv;
// const auto& blas_dhpmv = cblas_dhpmv;
// const auto& blas_dsymv = cblas_dsymv;
// const auto& blas_dsbmv = cblas_dsbmv;
// const auto& blas_dspmv = cblas_dspmv;
// const auto& blas_dtrmv = cblas_dtrmv;
// const auto& blas_dtbmv = cblas_dtbmv;
// const auto& blas_dtpmv = cblas_dtpmv;
// const auto& blas_dtrsv = cblas_dtrsv;
// const auto& blas_dtbsv = cblas_dtbsv;
// const auto& blas_dtpsv = cblas_dtpsv;

// const auto& blas_dger  = cblas_dger;
// const auto& blas_dgeru = cblas_dgeru;
// const auto& blas_dgerc = cblas_dgerc;
// const auto& blas_dher  = cblas_dher;
// const auto& blas_dhpr  = cblas_dhpr;
// const auto& blas_dher2 = cblas_dher2;
// const auto& blas_dhpr2 = cblas_dhpr2;
// const auto& blas_dsyr  = cblas_dsyr;
// const auto& blas_dspr  = cblas_dspr;
// const auto& blas_dsyr2 = cblas_dsyr2;
// const auto& blas_dspr2 = cblas_dspr2;

/*
 * BLAS Level 3 Routines (matrix-matrix)
 */

const auto& blas_dgemm = cblas_dgemm;
// const auto& blas_dsymm = cblas_dsymm;
// const auto& blas_dhemm = cblas_dhemm;
// const auto& blas_dsyrk = cblas_dsyrk;
// const auto& blas_dherk = cblas_dherk;
// const auto& blas_dsyr2k= cblas_dsyr2k;
// const auto& blas_dher2k= cblas_dher2k;
// const auto& blas_dtrmm = cblas_dtrmm;
// const auto& blas_dtrsm = cblas_dtrsm;
